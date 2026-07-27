Development guide
=================

Repository structure
--------------------

The Python package is organized by command category:

.. code-block:: text

   src/
     mdtbx/
       cli.py
       utils/
       build/
       trajectory/
       analysis/
       cv/
       agent/
   tests/
   example/
   pymol-plugins/

``mdtbx.cli`` scans the command category packages at startup. A module becomes a
subcommand when it exposes ``add_subcmd``. Helper modules without that function
are skipped.

Add a command
-------------

Place the module in the appropriate category and implement two functions:

.. code-block:: python

   def add_subcmd(subparsers):
       parser = subparsers.add_parser(
           "example",
           help="Describe one focused operation",
       )
       parser.add_argument("--input", required=True)
       parser.set_defaults(func=run)


   def run(args):
       ...

Keep the command focused on one operation. Reuse parsers and process helpers
from ``mdtbx.utils``. Avoid changing ``mdtbx.cli`` because registration is
automatic.

Add tests for parsing, normal behavior, failure behavior, and external command
construction. External tools should be mocked where practical.

Development commands
--------------------

.. code-block:: console

   $ pixi install
   $ pixi run test
   $ pixi run test-fast
   $ pixi run r
   $ pre-commit run --all-files

The full test suite can require packages from the locked environment but
should not require an actual production simulation.

Compatibility
-------------

The project currently fixes Python to version 3.10 and supports Linux x86-64,
macOS Apple silicon, and macOS Intel. Preserve command-line compatibility when
changing option names or defaults. Update the workflow guides, command tests,
and example scripts together when behavior changes.

Documentation coverage
----------------------

The command-reference test compares registered parser commands with
``:start_command:`` entries under ``docs/reference``. A new command is not
complete until it appears in the appropriate reference page and its user-facing
workflow documentation is updated when applicable.
