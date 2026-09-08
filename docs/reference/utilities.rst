Utility commands
================

cmd
---

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: cmd
   :groups:

convert
-------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: convert
   :groups:

mod_mdp
-------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: mod_mdp
   :groups:

partial_tempering
-----------------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: partial_tempering
   :groups:

rmfile
------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: rmfile
   :groups:

run_abfe
--------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: run_abfe
   :groups:

run_fep
-------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: run_fep
   :groups:

equilibrate_fep
---------------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: equilibrate_fep
   :groups:

shell_hook
----------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: shell_hook
   :groups:

skill
-----

Print the bundled English agent guide as Markdown to standard output. The guide
explains system preparation, simulation, analysis, and recovery, and directs
agents to follow the maintained workflows under ``example/``.

.. code-block:: console

   $ pixi run mdtbx skill
   $ pixi run mdtbx skill > MDTBX_SKILL.md

The guide is included in the installed package and can be displayed from any
working directory. Paths to examples in the guide refer to a source checkout.

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: skill
   :groups:

show_mdtraj
-----------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: show_mdtraj
   :groups:

show_npy
--------

.. autoprogram:: mdtbx.cli:create_parser()
   :prog: mdtbx
   :start_command: show_npy
   :groups:
