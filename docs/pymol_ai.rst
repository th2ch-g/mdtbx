PyMOL AI assistant
===================

After running ``pixi run pymolrc``, the PyMOL plugin exposes ``claude`` and
``codex`` commands. Each command captures a temporary screenshot and compact
scene metadata, requests PyMOL code, validates it, and applies it to the live
session.

.. code-block:: text

   PyMOL> claude color the ligand red
   PyMOL> codex show it as sticks and zoom to it
   PyMOL> ai_status
   PyMOL> ai_history
   PyMOL> ai_clear

Conversation behavior
---------------------

Completed turns are stored only in process memory. Claude and Codex share the
same history, so a request sent to one backend can refer to a result produced
by the other. The history is limited to 10 turns and 24 KiB and disappears when
PyMOL exits. Screenshots and scene metadata are not retained in history.

All requests use one worker queue and execute in submission order. The default
is asynchronous. Pass ``0`` as the second PyMOL argument to wait synchronously:

.. code-block:: text

   PyMOL> claude color chain A blue, 0

``ai_status`` shows queued, running, retrying, executing, completed, and failed
jobs. ``ai_history`` prints retained completed turns; an integer argument limits
the output. ``ai_clear`` clears conversation history without changing the
scene or deleting job status.

Safety boundary
---------------

Generated code is still auto-applied. Python is restricted to public
``cmd.*`` calls and a small built-in set. Imports, private attributes, shell
execution, script loading, Python-evaluation commands, and unknown PyMOL
commands are rejected. A failed generated response is returned to the backend
for correction, up to the configured attempt limit.

Claude runs without tools or interactive permission prompts. Codex runs with
an ephemeral session, a read-only sandbox, and approvals disabled. These
controls constrain the AI CLI; they do not replace review of the resulting
live molecular visualization.
