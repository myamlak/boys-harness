# Registration check

`check_registration.py` checks the documented-claims enumeration against the
registration table in the boys accuracy gate, so the two cannot drift.

## What it reads

- `--boys <checkout>`: the gate source at
  `<checkout>/tests/boys_accuracy_gate.cpp`.
- `--tsv <path>`: the enumeration. The default is
  `accuracy-claim-audit/records/documented-claims-v5.tsv`, resolved from this
  directory.

From the gate source it parses, under the preprocessor configuration the
published claims-book numbering was made in (`BOYS_GATE_FP16` defined,
`BOYS_GATE_EXTERNAL_BOOST` not defined):

- the claims book's `add("row", ...)` call sites, in printed order — the
  enumeration's G indices;
- `kRegistrations[]`, one entry per documented claim with the row it names;
- `kUnregisteredRows[]`, the rows the gate marks as carrying no registered
  claim;
- the book's row identifiers and their order.

From the enumeration it reads the rows whose status is `claims-book` and their
`registered_by` G indices.

## The checks

(a) The claim ids with status `claims-book` in the enumeration and the claim
ids in `kRegistrations[]` are the same set.

(b) For every claims-book claim:

- every `registered_by` G index is a row of the claims book;
- the row the gate registers for the claim is one of the rows the claim cites;
- every row the claim cites is registered to the claim, registered to nobody,
  or registered to another claim that also cites it.

The enumeration's numbering is frozen at the book it was made from. A
registration to a row appended after the enumeration's highest cited index
cannot be cited; it is accepted only when every row the claim cites is marked
in `kUnregisteredRows` or is cited by another claim registered to it, and is
printed as a `notice:` with the claim, the cited rows and the appended row.
`--strict` counts each such registration as a mismatch. Once the enumeration
cites the appended row, the ordinary citation rule applies to it.

(c) Every row named by a registration is a row of the claims book.

The script prints the G-index mapping it used, with the source line of each
`add()` call.

## The command

```sh
python accuracy-claim-audit/tools/check_registration.py --boys <boys-checkout>
```

`--tsv <path>` checks another copy of the enumeration. `--strict` treats a
registration to a row appended after the enumeration's reach as a mismatch.

Exit status 0 when the two agree, 2 on any mismatch or parse failure. A
checkout whose gate has no registration table is refused with status 2 and a
message naming the file; so is a source or enumeration that cannot be parsed.

Only the Python standard library is used.
