The Spectral Sequences Project
==============================

This is a monorepo containing various projects:

1. `ext`
A general library to work with Ext over an Fp algebra. More generally, it
allows us to work compute in the derived category of said algebra. The primary
purpose is to compute the classical Adams E2 page by computing Ext over the
Steenrod algebra.

**To get started computing Ext over the Steenrod algebra, read the Quickstart [here](ext/README.md).**

2. `web_ext`
Web interfaces to `ext`. There are two subprojects at the moment:

 - `sseq_gui`: A GUI to work with the Adams spectral sequence. Given a
   Steenrod module, this computes its Ext and displays the associated Adams
   spectral sequence. The user can then interactively input differentials and
   the program can propagate differentials via the Leibniz rule.

   This can be tried out at https://spectralsequences.github.io/sseq/ which
   does not require installation.

 - `steenrod_calculator`: This is a simple user interface to compute sums and
   products in the Steenrod algebra and express the result in either the Adem
   or Milnor basis.

   This is available at
   https://spectralsequences.github.io/sseq/calculator/ .

 - `module_builder`: An interactive editor for finite dimensional modules over
   the Steenrod algebra. Add cells, draw the action of the algebra generators
   between them, and have the Adem relations checked as you edit. Modules save
   in the same format as `ext/steenrod_modules`, and can be handed straight to
   `sseq_gui`.

   This is available at
   https://spectralsequences.github.io/sseq/module_builder/ .

3. `python_ext`
WIP python bindings for the `ext` library.

4. `chart`
A general spectral sequence web interface, with a python-based repl for
programmatic interaction.
