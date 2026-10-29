# Great Pumpkin: a small Lean 4 experiment

The [mathematical note](../../docs/expeditions/great-pumpkin-001-mathematics.md)
states the conventions, proofs and their physical limits. The project needs
Lean **4.24.0**, with its bundled standard library; no Mathlib or other Lake
dependency. Do not substitute an unspecified latest toolchain for a reproducible
check. Install Lean using its [official instructions](https://lean-lang.org/install/),
or use the [pinned release](https://github.com/leanprover/lean4/releases/tag/v4.24.0).

From this directory:

```text
lake build
lake env lean --trust=0 GreatPumpkin.lean
```

The second invocation elaborates the source in a fresh process and asks Lean
to type-check imported modules too. It does not rely on our previously built
project `.olean`. This is an independent kernel invocation, not an independent
implementation of the Lean kernel. Check that all nine `#print axioms` reports
list only `propext`, `Classical.choice`, and `Quot.sound`.

To record a versioned, source-hashed proof-check report from repository root:

```text
python experiments/great_pumpkin_001/check_formal.py --lean PATH_TO_LEAN_EXECUTABLE
```

The expedition used the Windows x64 official release archive with SHA-256
`c2fcf1cba1089e526e184de2a6f6763a6e15596a2991f096dfacd4a6e3ce2cdd`.
This is a recorded download hash, not an independently authenticated signature.
The compiler reported commit `797c613eb9b6d4ec95db23e3e00af9ac6657f24b`.
Tool binaries and build products are ignored, not vendored into this repository.

The theorems quantify over a type with standard ordered-field instances; their
laws are explicit assumptions and are instantiated by `Rat` in the file.
They apply algebraically to an ordered real field as well, but the project
does not construct the reals or formalize the physical logarithm mapping.
Completeness, quantum field theory and empirical input accuracy are not
smuggled in as axioms.
