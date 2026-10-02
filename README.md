# Qwen and protein-model alignment dashboard

Static browser-side report for a fixed panel of 2,048 proteins, five Qwen checkpoints, and nine protein-model checkpoints. The site has no sign-in, analytics, external/CDN dependencies, or analysis backend. It downloads its bundled assets from the same static host and runs the report in the browser. Imported dashboard state remains in the browser.

## Provenance

- Hosting package date: 2026-10-02.
- Analysis/report snapshot: alignment review 2026-10-01.
- Original standalone HTML SHA-256: `c0704cdc09efae87c45ac1c0a8d6d44f28bf227066bae04c7b1dd47e1aef3c05`.
- Original `dashboard_data.json` SHA-256: `c75cc4ddf34ccd285de6cbdf9d9becb1c11bc9b5fac31031a453e73d1b8b4f00`.
- Decompressing `assets/dashboard_data.json.gz` restores the original JSON bytes exactly. No values, protein IDs, cohort rules, metrics, or plot layouts were changed for hosting.
- Hosting-only changes: external local assets, compressed-data loader/loading and error display, and network/offline wording appropriate for a website.

## Loading

`index.html` loads local CSS and the loader. The loader downloads the gzip data and Plotly, verifies the JSON checksum when browser cryptography is available, then starts `dashboard_ui.js`, `dashboard_cka.js`, and `dashboard_descriptive.js` in that order. Native `DecompressionStream` is preferred; bundled pako handles browsers without it. This hosted version must be opened over HTTP(S); the original standalone report remains the version for opening directly from disk.

## Third-party notices

- Plotly.js 3.1.0 is the exact distribution embedded in the frozen report. Its MIT license notice is retained in `assets/plotly.min.js`.
- pako 1.0.11 is vendored from the existing local runtime solely for gzip decompression. Its MIT/Zlib attribution is in `assets/pako-LICENSE.txt` and the distribution header.
