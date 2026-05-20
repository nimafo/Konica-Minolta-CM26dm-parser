# Spectral material parser
## Konica Minolta CM26dm parser
A lightweight Python class for processing reflectance data from:
1. The **Konica Minolta CM-26dG** spectrophotometer.
2. The [SpectralDB](https://github.com/C38C/SpectralDB).
   
---

## Features

- Parses CM-26dG CSV exports (SCI + SCE measurement pairs)
- Converts CIE L\*a\*b\* (D65) → XYZ → sRGB via the [`colour`](https://colour.science) library
- Computes visible reflectance
- Estimates surface specularity from the SCI–SCE difference

---

## Requirements

```
pandas
numpy
colour-science
```

Install with:

```bash
pip install pandas numpy colour-science
```

---

## Usage

```python
from GT import GT

# Load material #3 from a CM-26dG export
mat = GT("20241004_image2mat_640.csv", n=3)

print(mat.Lab_SCI)          # CIE L*a*b* (SCI, D65)
print(mat.RGB_SCI)          # sRGB [0–1]
print(mat.reflectance_SCI)  # Visible reflectance [0–1]
print(mat.specularity)      # SCI − SCE reflectance difference
print(mat.spectral_SCI)     # Full spectral curve (360–740 nm)
```

For SpectralDB:

```python
from GT import GT_fromSpectralDB

mat = GT_fromSpectralDB("spectral_db.csv", n=1)
```

---

## CSV Format

The input CSV must follow the CM-26dG export structure: **two rows per material** (one `SCI`, one `SCE`), with columns `L*(D65)`, `a*(D65)`, `b*(D65)`, `Target Name`, and wavelength columns `360nm` through `740nm`.

---

## Attributes

| Attribute | Description |
|---|---|
| `material_count` | Total number of materials in the file |
| `spectral_SCI/SCE` | Spectral reflectance Series, 360–740 nm |
| `Lab_SCI/SCE` | CIE L\*a\*b\* values under D65 |
| `RGB_SCI/SCE` | sRGB values |
| `reflectance_SCI/SCE` | visible reflectance |
| `specularity` | Surface specularity |
