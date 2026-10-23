# MAX identities with two hidden ReLU layers

Certificates and verification code for:

**Zhimao Wang and Amitabh Basu. [Representing MAX functions using two-hidden-layer ReLU networks](https://arxiv.org/abs/2608.25221). arXiv:2608.25221, 2026.**

The certificates express $`2\mathrm{MAX}_N(x)=2\max\{x_1,\ldots,x_N\}`$ as rational linear combinations of nested maximum expressions, for $N=5,6,7,8$. Each expression has a two-hidden-layer ReLU realization. Verification uses the authors' rank-counting method.

## Run the verifier

Requires Python 3; no external packages. From the repository root:

```bash
python verify_MAX_identity.py certificates/identity_MAX_7_k3.txt
```

Each run checks one file. Omit the filename for an interactive prompt. Files in `negative_examples/` have deliberately altered coefficients and are expected to fail.

## Certificate notation

Each file begins with a target such as `2 MAX_7 =`, followed by terms of the form `coefficient * T[left | right]`.

A word `ab` denotes $\max(x_a,x_b)$, and `aa` denotes $x_a$. For example, `T[aa ab | cd ef]` sums expressions of the form

$$
\max\{x_a+\max(x_a,x_b),\;\max(x_c,x_d)+\max(x_e,x_f)\}
$$

over **all injective assignments of the distinct pattern letters to the N coordinates**, including symmetry multiplicities. Each side has $k$ pairwise maxima, with the same $k$ throughout the certificate.

## What the verifier checks

By symmetry, verification takes place on $x_1\geq\cdots\geq x_N$. Outer maxima with a side that dominates throughout this cone contribute to the fixed linear vector; the remaining maxima are grouped by their pairs of coefficient vectors as ambiguous classes.

Using integer counts and exact `Fraction` arithmetic, the verifier checks that **both the fixed residual vector and all ambiguous coefficients vanish** after subtracting the target.

`VALID IDENTITY` proves the stated identity. `INVALID IDENTITY` means this cancellation check failed; it need not disprove a general input identity, since further merging of distinct ambiguous classes with the same direction is not performed.

## Citation

```bibtex
@misc{wang2026representing,
  title={Representing MAX functions using two-hidden-layer ReLU networks},
  author={Zhimao Wang and Amitabh Basu},
  year={2026},
  eprint={2608.25221},
  archivePrefix={arXiv},
  url={https://arxiv.org/abs/2608.25221}
}
```
