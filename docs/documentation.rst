Documentation and translation
=============================

Source layout
-------------

English reStructuredText files under ``docs/`` are the canonical source.
Japanese translations are stored as gettext catalogs under
``docs/locale/ja/LC_MESSAGES/``. Generated ``.pot`` and ``.mo`` files are not
committed.

Build locally
-------------

Install and build the documentation environment:

.. code-block:: console

   $ pixi install -e docs
   $ pixi run -e docs docs

The English site is written to ``docs/_build/html/`` and the Japanese site to
``docs/_build/html/ja/``.

Update Japanese catalogs
------------------------

After editing English source:

.. code-block:: console

   $ pixi run -e docs docs-update-ja

Translate new or changed ``msgid`` entries in the corresponding ``.po`` files.
Do not translate commands, file names, option names, paths, or code samples.
Preserve reStructuredText roles and cross-reference targets exactly.

Validate changes
----------------

.. code-block:: console

   $ pixi run -e docs docs
   $ pixi run -e docs docs-linkcheck
   $ pixi run test

Both HTML builders treat warnings as errors. Check the English and Japanese
navigation, search, language links, command options, and code blocks before
merging.
