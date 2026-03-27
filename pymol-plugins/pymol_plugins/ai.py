"""AI assistant command for PyMOL.

Usage:
    claude <natural language instruction>
    codex <natural language instruction>

Examples:
    claude resid 10を赤色にして
    codex タンパク質を透明度0.5で表示して
"""

from __future__ import annotations

import base64
import json
import re
import subprocess
import tempfile
import threading
import traceback
from itertools import count
from pathlib import Path

from pymol import cmd

SUPPORTED_AI_TYPES = {"claude", "codex"}
AI_TIMEOUT_SEC = 180
AI_JOB_COUNTER = count(1)
AI_JOBS: dict[int, dict[str, object]] = {}
AI_JOBS_LOCK = threading.Lock()


def _capture_scene_png(image_path: Path) -> None:
    """Capture the current viewport using the PyMOL `png` command."""
    viewport = tuple(cmd.get_viewport())
    try:
        cmd.refresh()
        cmd.png(str(image_path), ray=0, quiet=1)
    finally:
        if len(viewport) == 2:
            cmd.viewport(int(viewport[0]), int(viewport[1]))
        cmd.refresh()

    if not image_path.exists():
        raise RuntimeError(f"failed to save screenshot to {image_path}")


def _get_scene_context() -> str:
    """Collect compact PyMOL scene metadata for the AI prompt."""
    objects = cmd.get_object_list("all")
    names = cmd.get_names("all")
    enabled_objects = cmd.get_names("objects", enabled_only=1)

    lines = [
        f"Loaded objects: {objects}",
        f"Enabled objects: {enabled_objects}",
        f"All named selections/objects: {names}",
    ]

    for obj in objects[:3]:
        residues: list[tuple[str, str, str]] = []
        cmd.iterate_state(
            1,
            f"({obj}) and name CA",
            "residues.append((resi, resn, chain))",
            space={"residues": residues},
        )
        if not residues:
            continue

        summary = ", ".join(
            f"{chain or '-'}:{resn}{resi}" for resi, resn, chain in residues[:10]
        )
        if len(residues) > 10:
            summary += f" ... ({len(residues)} residues total)"
        lines.append(f"{obj} residues: {summary}")

    return "\n".join(lines)


def _build_prompt(instruction: str, scene_context: str, image_path: Path) -> str:
    """Build the prompt sent to the local AI CLI."""
    return f"""You are a PyMOL expert assistant.

The user is controlling a live PyMOL session. A PNG screenshot of the current viewport is attached.
The screenshot file name is `{image_path.name}`.

Current PyMOL scene metadata:
{scene_context}

User request:
{instruction}

Return only executable PyMOL commands or Python code for PyMOL.

Rules:
- No explanation, no prose, no markdown outside code fences.
- Prefer PyMOL commands when simple enough.
- If you return Python, use a single ```python fenced block and call `cmd.*`.
- If you return PyMOL commands, use a single ```pymol fenced block.
- Keep the output minimal and directly executable in the current session.
- Do not use shell commands.
"""


def _normalize_bool_arg(value: object, default: bool = True) -> bool:
    if value is None:
        return default

    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    return default


def _extract_code(response: str) -> list[tuple[str, str]]:
    """Extract executable code blocks from the AI response."""
    blocks: list[tuple[str, str]] = []
    pattern = re.compile(r"```(pymol|python)\s*\n(.*?)```", re.DOTALL | re.IGNORECASE)

    for match in pattern.finditer(response):
        lang = match.group(1).lower()
        code = match.group(2).strip()
        if code:
            blocks.append((lang, code))

    if blocks:
        return blocks

    stripped = response.strip()
    if stripped and _looks_like_raw_pymol_commands(stripped):
        return [("pymol", stripped)]
    return []


def _run_subprocess(
    args: list[str], stdin: str | None = None
) -> subprocess.CompletedProcess:
    return subprocess.run(
        args,
        input=stdin,
        capture_output=True,
        text=True,
        timeout=AI_TIMEOUT_SEC,
        check=False,
    )


def _collect_text_from_message_content(content: object) -> str:
    if isinstance(content, str):
        return content

    if not isinstance(content, list):
        return ""

    texts: list[str] = []
    for block in content:
        if not isinstance(block, dict):
            continue
        if block.get("type") == "text" and isinstance(block.get("text"), str):
            texts.append(block["text"])
    return "".join(texts)


def _parse_claude_stream_json(stdout: str) -> str:
    """Convert Claude stream-json output into plain assistant text."""
    final_texts: list[str] = []
    streamed_text_parts: list[str] = []

    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue

        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            final_texts.append(line)
            continue

        event_type = event.get("type")

        if event_type == "result":
            result_text = event.get("result")
            if isinstance(result_text, str) and result_text.strip():
                final_texts.append(result_text)

            message = event.get("message")
            if isinstance(message, dict):
                text = _collect_text_from_message_content(message.get("content"))
                if text:
                    final_texts.append(text)
            continue

        if event_type == "assistant":
            message = event.get("message")
            if isinstance(message, dict):
                text = _collect_text_from_message_content(message.get("content"))
                if text:
                    final_texts.append(text)
            continue

        if event_type == "stream_event":
            raw_event = event.get("event")
            if not isinstance(raw_event, dict):
                continue

            if raw_event.get("type") != "content_block_delta":
                continue

            delta = raw_event.get("delta")
            if not isinstance(delta, dict):
                continue

            if delta.get("type") == "text_delta" and isinstance(delta.get("text"), str):
                streamed_text_parts.append(delta["text"])

    if final_texts:
        return "\n".join(text for text in final_texts if text.strip()).strip()
    return "".join(streamed_text_parts).strip()


def _has_executable_content(lang: str, code: str) -> bool:
    if lang == "python":
        if not code.strip():
            return False
        try:
            compile(code, "<pymol-ai>", "exec")
        except SyntaxError:
            return False
        return True

    return any(
        line.strip() and not line.strip().startswith("#") for line in code.splitlines()
    )


def _looks_like_raw_pymol_commands(text: str) -> bool:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return False

    allowed_prefixes = (
        "@",
        "/",
        "run ",
        "python",
        "python end",
        "embed",
        "skip",
    )
    common_commands = {
        "align",
        "bg_color",
        "cartoon",
        "center",
        "color",
        "delete",
        "disable",
        "dist",
        "distance",
        "extract",
        "fetch",
        "hide",
        "label",
        "load",
        "orient",
        "png",
        "ray",
        "remove",
        "rebuild",
        "refresh",
        "reinitialize",
        "select",
        "set",
        "show",
        "spectrum",
        "super",
        "turn",
        "util.cbag",
        "util.cbaw",
        "util.cbao",
        "zoom",
    }

    for line in lines:
        if line.startswith("#"):
            continue
        if line.startswith(allowed_prefixes):
            continue

        first = line.split(maxsplit=1)[0].rstrip(",").lower()
        if first in common_commands:
            continue
        return False

    return True


def _run_claude(prompt: str, image_path: Path) -> str:
    """Run Claude Code in stream-json mode to attach an image."""
    image_b64 = base64.b64encode(image_path.read_bytes()).decode("ascii")
    payload = {
        "type": "user",
        "message": {
            "role": "user",
            "content": [
                {"type": "text", "text": prompt},
                {
                    "type": "image",
                    "source": {
                        "type": "base64",
                        "media_type": "image/png",
                        "data": image_b64,
                    },
                },
            ],
        },
    }

    result = _run_subprocess(
        [
            "claude",
            "--print",
            "--verbose",
            "--input-format",
            "stream-json",
            "--no-session-persistence",
            "--output-format",
            "stream-json",
        ],
        stdin=json.dumps(payload) + "\n",
    )
    if result.returncode != 0:
        stderr = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"claude command failed (exit {result.returncode}):\n{stderr}"
        )
    return _parse_claude_stream_json(result.stdout)


def _run_codex(prompt: str, image_path: Path) -> str:
    """Run Codex CLI with an attached screenshot and capture the last message."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "codex-last-message.txt"
        result = _run_subprocess(
            [
                "codex",
                "exec",
                "--skip-git-repo-check",
                "-c",
                "effort=high",
                "--image",
                str(image_path),
                "--output-last-message",
                str(output_path),
                prompt,
            ]
        )
        if result.returncode != 0:
            stderr = result.stderr.strip() or result.stdout.strip()
            raise RuntimeError(
                f"codex command failed (exit {result.returncode}):\n{stderr}"
            )

        if output_path.exists():
            return output_path.read_text().strip()

    return result.stdout


def _register_ai_job(job_id: int, ai_type: str, instruction: str) -> None:
    with AI_JOBS_LOCK:
        AI_JOBS[job_id] = {
            "id": job_id,
            "type": ai_type,
            "instruction": instruction,
            "status": "running",
        }


def _update_ai_job(job_id: int, **updates: object) -> None:
    with AI_JOBS_LOCK:
        job = AI_JOBS.get(job_id)
        if job is None:
            return
        job.update(updates)


def _snapshot_ai_jobs() -> list[dict[str, object]]:
    with AI_JOBS_LOCK:
        return [dict(job) for job in AI_JOBS.values()]


def _execute_blocks(blocks: list[tuple[str, str]]) -> None:
    """Execute extracted code blocks inside the active PyMOL session."""
    for lang, code in blocks:
        if not _has_executable_content(lang, code):
            continue

        print(f" [ai] Executing ({lang}):\n{code}\n")
        if lang == "python":
            exec(code, {"cmd": cmd, "__builtins__": __builtins__})  # noqa: S102
            continue

        for line in code.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                cmd.do(line)


def _run_ai_job(
    job_id: int,
    ai_type: str,
    instruction: str,
    prompt: str,
    image_path: str,
) -> None:
    try:
        response = (
            _run_claude(prompt, Path(image_path))
            if ai_type == "claude"
            else _run_codex(prompt, Path(image_path))
        )
        response = response.strip()
        _update_ai_job(job_id, status="completed", response=response)
        print(f" [ai:{job_id}] Response received ({len(response)} chars).")

        blocks = [
            (lang, code)
            for lang, code in _extract_code(response)
            if _has_executable_content(lang, code)
        ]
        if not blocks:
            _update_ai_job(job_id, status="no_code")
            print(f" [ai:{job_id}] No executable code found in response.")
            print(f" [ai:{job_id}] AI response:\n{response}")
            return

        _update_ai_job(job_id, status="executing", blocks=len(blocks))
        print(f" [ai:{job_id}] Executing {len(blocks)} block(s).")
        _execute_blocks(blocks)
        _update_ai_job(job_id, status="done")
        print(f" [ai:{job_id}] Done.")
    except FileNotFoundError:
        _update_ai_job(job_id, status="error", error=f"'{ai_type}' command not found")
        print(
            f" [ai:{job_id}] Error: '{ai_type}' command not found. "
            "Make sure it is installed and available in PATH."
        )
    except subprocess.TimeoutExpired:
        _update_ai_job(job_id, status="error", error=f"timeout ({AI_TIMEOUT_SEC}s)")
        print(f" [ai:{job_id}] Error: AI command timed out ({AI_TIMEOUT_SEC}s).")
    except RuntimeError as exc:
        _update_ai_job(job_id, status="error", error=str(exc))
        print(f" [ai:{job_id}] Error: {exc}")
    except Exception as exc:  # noqa: BLE001
        error = f"{type(exc).__name__}: {exc}"
        _update_ai_job(job_id, status="error", error=error)
        print(f" [ai:{job_id}] Error while applying response: {error}")
        print(traceback.format_exc().rstrip())
    finally:
        try:
            Path(image_path).unlink(missing_ok=True)
            Path(image_path).parent.rmdir()
        except OSError:
            pass


def _submit_ai_request(
    instruction: str,
    type: str = "claude",  # noqa: A002
    async_: str = "1",
) -> None:
    """Send a natural-language request to Claude/Codex and run the result."""
    ai_type = type.strip().lower()
    if ai_type not in SUPPORTED_AI_TYPES:
        supported = ", ".join(sorted(SUPPORTED_AI_TYPES))
        print(f" [ai] Unknown type '{type}'. Use one of: {supported}")
        return

    if not instruction or not instruction.strip():
        print(" [ai] Instruction is empty.")
        return

    use_async = _normalize_bool_arg(async_, default=True)
    image_dir = Path(tempfile.mkdtemp(prefix="pymol-ai-"))
    image_path = image_dir / "pymol_scene.png"
    _capture_scene_png(image_path)
    scene_context = _get_scene_context()
    prompt = _build_prompt(instruction.strip(), scene_context, image_path)

    job_id = next(AI_JOB_COUNTER)
    _register_ai_job(job_id, ai_type, instruction.strip())

    print(
        f" [ai:{job_id}] Using {ai_type}. "
        f"{'Queued asynchronously' if use_async else 'Running synchronously'}: {instruction!r}"
    )
    print(f" [ai:{job_id}] Screenshot saved to {image_path}")

    if use_async:
        worker = threading.Thread(
            target=_run_ai_job,
            args=(job_id, ai_type, instruction.strip(), prompt, str(image_path)),
            daemon=True,
            name=f"pymol-ai-{job_id}",
        )
        _update_ai_job(job_id, thread=worker)
        worker.start()
        print(
            f" [ai:{job_id}] Background job started. Use `ai_status` to check progress."
        )
        return

    _run_ai_job(job_id, ai_type, instruction.strip(), prompt, str(image_path))


def claude_cmd(instruction: str, async_: str = "1") -> None:
    """Run the local Claude CLI for a PyMOL request."""
    _submit_ai_request(instruction=instruction, type="claude", async_=async_)


def codex_cmd(instruction: str, async_: str = "1") -> None:
    """Run the local Codex CLI for a PyMOL request."""
    _submit_ai_request(instruction=instruction, type="codex", async_=async_)


def ai_status(job_id: str = "all") -> None:
    """Show status of AI jobs."""
    jobs = _snapshot_ai_jobs()
    if not jobs:
        print(" [ai] No jobs.")
        return

    if job_id.strip().lower() != "all":
        try:
            wanted = int(job_id)
        except ValueError:
            print(" [ai] job_id must be an integer or 'all'.")
            return
        jobs = [job for job in jobs if job.get("id") == wanted]
        if not jobs:
            print(f" [ai] No such job: {wanted}")
            return

    for job in sorted(jobs, key=lambda item: int(item["id"])):
        status = job.get("status", "unknown")
        ai_type = job.get("type", "?")
        instruction = job.get("instruction", "")
        print(f" [ai:{job['id']}] {status} ({ai_type}) {instruction}")
        if status == "error" and job.get("error"):
            print(f" [ai:{job['id']}] error: {job['error']}")


cmd.extend("claude", claude_cmd)
cmd.extend("codex", codex_cmd)
cmd.extend("ai_status", ai_status)
