# Kernel release identity across v5.3.1, v5.3.2 and v5.3.3

This record answers the question the paper's version-bound evidence raises: the
consumer session was captured at v5.3.1, the libint2 probe was built against
v5.3.2 headers, and the release is v5.3.3. Are the evaluated files the same
arithmetic at all three?

The check is mechanical. It compares every file under `include/boys/` and
`src/` at the three annotated tags, and for the one source file that differs it
compares a comment-stripped canonical form.

## Tags

| tag | commit |
| --- | --- |
| v5.3.1 | `2f83e40a6504368b7b2f13682e29b4f24ab2a40c` |
| v5.3.2 | `6739447d82765ae5f7f68351a99f418cd018dd9a` |
| v5.3.3 | `01bba7b0ed7e0dd387f2af087e992445222e5546` |

## Command

From a clone of `github.com/myamlak/boys`:

```sh
for t in v5.3.1 v5.3.2 v5.3.3; do
  for f in $(git ls-tree -r --name-only "$t" include/boys src); do
    printf '%s  %s  %s\n' "$t" "$(git show "$t:$f" | sha256sum | cut -d' ' -f1)" "$f"
  done
done
```

`identity.json` carries the result. The comment-stripped hash of
`src/boys_simd.cpp` was computed with `check_source_identity.py` in this
directory.

## Result

- 78 of the 80 files are byte-identical at all three tags.
- `include/boys/version.hpp` differs at each tag: it carries the version
  string, not arithmetic.
- `src/boys_simd.cpp` differs between v5.3.1/v5.3.2 and v5.3.3. The difference
  is comments only. Its comment-stripped canonical hash is
  `5f757ee418140e3cceeeb059fb2dc65ad9077815334be3501385a05ffd164ca8` at all
  three tags.
- No other file under `include/boys/` or `src/` differs.

The libint2 probe includes the public headers and the consumer session links
the library built from these sources. The evaluated arithmetic is therefore the
release arithmetic. The gate at v5.3.3 is separately run on its committed grid,
and the exhaustive half sweep has its own v5.3.3 re-run in
`../half-sweep/v533-rerun/`.

## Scope

This is a source-identity check, not a bit-for-bit binary comparison. A
different compiler or flag set can change code generation, so the record is
stated for the build settings the paper's measurements use. It closes the
version-carrying question for the source; it does not replace the gate run or
the half-sweep re-run.
