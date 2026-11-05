# Module 03 synthetic classification test

Run the synthetic end-to-end test from the repository root:

```bash
python scripts/03_pan_TE_catalogue/tests/test_classification_release.py
```

The test creates temporary synthetic records and exercises the evidence-classification and rescue stages. It is a software regression test; it does not use the study's production TE libraries and does not validate biological classification accuracy on the original assemblies.

GitHub Actions runs this test when this directory or the classification scripts change, including on the public-release commit.
