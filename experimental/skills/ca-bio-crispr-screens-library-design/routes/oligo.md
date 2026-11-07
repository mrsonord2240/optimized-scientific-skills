# Format an existing spacer
~~~bash
python scripts/build_oligo.py "$SPACER" --subpool 1
~~~
Set SPACER to an existing reviewed 20-base sequence. The source helper concatenates its built-in subpool primer and scaffold stub; it does not validate sequence length or alphabet, and its --vector option does not change the result. Verify sequence, orientation, intended architecture and synthesis requirements independently. See [source format caveats](../references/oligo-design.md).
Done when the output is checked against the intended input and architecture.
