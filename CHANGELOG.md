# Changelog

## 0.5.0

- Add optional `compile_qram`: shared phase-gradient variable rotations and explicit four-word read-only Boolean QRAM.
- Preserve full controlled complex phase, clean workspace, public interfaces, and compressed calls/repeats.
- Reject write-bearing/open programs and reserved namespace conflicts; count gradient preparation and cleanup.
- Add interpreter and both native PySparQ adapter equivalence checks; document resource and timing tradeoffs in both languages.
