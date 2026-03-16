# chemkin2suth

Convert Chemkin transport files into OpenFOAM Sutherland transport coefficients that can be read by **ChemkinToFoam**

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://xxxx
cd chemkin2suth
pip install -e .
```

Required Python version:

Python >= 3.9

Dependencies will be installed automatically:

- numpy
- scipy

---

## Command Line Usage

Convert a Chemkin transport file:

```bash
chemkin2suth transport.dat thermo.dat transportProperties
```

### Arguments

| Argument | Description |
|------|------|
| transport.dat | Chemkin transport input file |
| thermo.dat | Chemkin thermo input file |
| transportProperties | Output OpenFOAM Sutherland transport dictionary |

Example:

```bash
chemkin2suth transport.dat thermo.dat transportDict
```

The output file will contain entries like:

```
AR
{
    transport
    {
        As              1.716e-05;
        Ts              111;
    }
}
```

---

## Python API

The converter can also be used inside Python:

```python
from chemkin2suth import convert_transport_file

convert_transport_file(
    "transport.dat",
    "thermo.dat",
    "transportDict"
)
```

---

## Method

For each species:

1. Read Chemkin transport parameters
2. Compute viscosity over a temperature range
3. Fit the Sutherland law

μ(T) = As * T^(3/2) / (T + Ts)

4. Output fitted coefficients for OpenFOAM.

---

## Project Structure

```
chemkin2suth/
├── pyproject.toml
├── README.md
├── src/
│   └── chemkin2suth/
│       ├── __init__.py
│       ├── cli.py
│       └── converter.py
```

---

## License

MIT License

---

## Author

Zijian Sun  
Princeton University