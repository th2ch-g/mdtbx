#!/usr/bin/env python3
"""Claude Code approval boundary for mdtbx Agent operations."""

from __future__ import annotations

import json
import re
import sys

MUTATING_SCHEDULER = re.compile(
    r"(?<![\w/.-])(sbatch|qsub|pjsub|scancel|qdel|pjdel)(?![\w.-])"
)
APPROVAL_COMMAND = re.compile(
    r"(?<![\w/.-])mdtbx\s+(agent_run|agent_profile_save)(?![\w.-])"
)


def decision(value: str, reason: str) -> dict:
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": value,
            "permissionDecisionReason": reason,
        }
    }


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, OSError):
        return 0
    command = payload.get("tool_input", {}).get("command", "")
    if not isinstance(command, str):
        return 0
    if MUTATING_SCHEDULER.search(command):
        print(
            json.dumps(
                decision(
                    "deny",
                    "Use an immutable mdtbx agent_plan and approved agent_run.",
                )
            )
        )
        return 0
    match = APPROVAL_COMMAND.search(command)
    if match:
        print(
            json.dumps(
                decision(
                    "ask",
                    f"Approve the exact mdtbx {match.group(1)} hash in this command.",
                )
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
