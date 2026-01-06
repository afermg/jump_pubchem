import marimo

__generated_with = "0.11.31-dev3"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import duckdb
    import pandas as pd
    import pyarrow
    from pyodide.http import pyfetch
    return duckdb, mo, pd, pyarrow, pyfetch


@app.cell
def _(pd):
    pq = pd.read_parquet('https://zenodo.org/api/records/15104315/files/jump_pubchem.parquet/content')
    return (pq,)


@app.cell
def _(mo):
    mapper = mo.sql(
        f"""
        SELECT *,(((from_hex(fingerprint)::BITSTRING)::VARCHAR)[33:-8]) AS bs FROM pq;
        """,
        output=False
    )
    return (mapper,)


@app.cell
def _(mo):
    sb = mo.ui.form(mo.ui.text(value=""))
    name_fields = {"Common Name":"name","InChIKey":"inchikey","SMILES":"smiles"}
    field = mo.ui.dropdown(options=name_fields,value=list(name_fields)[0])
    return field, name_fields, sb


@app.cell
def _(mo):
    mo.md(
        r"""
        # Find your compound in JUMP

        This tool uses Pubchem fingerprints to associate any Pubchem-available compound to its closest JUMP analog.

        Instructions:

        1. Submit your compound and identifier type (e.g., 'aspirin' and 'Common Name', 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N' and 'InChIKey')
        2. Copy the top choice (either InChiKey or Metadata_JCP2022) and use it on [broad.io/compound](http://broad.io/compound).
        3. If your compound is found on PubChem, you will be shown the top matches in JUMP
        """
    )
    return


@app.cell
def _(field, mo, sb):
    mo.hstack((mo.vstack((mo.md("Identifier"),sb)), mo.vstack((mo.md("Identifier type"), field))), justify="start")
    return


@app.cell
async def _(field, mo, pyfetch, sb):
    if sb.value:
        sim_col = "similarity"
        # result = get_pubchem_id(sb.value, field.value)
        result = await pyfetch(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{field.value}/{sb.value}/json")
        fmt = await result.json()
        if "Fault" in fmt:
            output = mo.md("Error: ID not found")
        else:
            options = [{x["urn"]["label"]:list(x["value"].values())[0] for x in hit["props"] if (x["urn"]["label"]=="InChIKey") or (x["urn"]["label"] == "IUPAC Name" and x["urn"]["name"]=="Preferred") or (x["urn"]["label"]=="Fingerprint") } for hit in fmt["PC_Compounds"]]
            # Convert to bitstring
            fp_hex = options[0]["Fingerprint"]
            bitstring = (bin(int(fp_hex, 16))).zfill(len(fp_hex)*4)[32:-7]

            # Query against JUMP->PubChem mapper
            output = mo.sql(
                f"SELECT Metadata_JCP2022,InChIKey,BIT_COUNT(bs::BITSTRING & '{bitstring}') / BIT_COUNT(bs::BITSTRING | '{bitstring}') AS {sim_col} FROM mapper ORDER BY {sim_col} DESC LIMIT 3"
            )
        output
    return bitstring, fmt, fp_hex, options, output, result, sim_col


if __name__ == "__main__":
    app.run()
