import marimo

__generated_with = "0.11.23"
app = marimo.App(width="medium")


@app.cell
def _():
    import json

    import duckdb
    import marimo as mo
    import requests

    def get_pubchem_id(common_name, identifier):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{identifier}/{common_name}/json"
        response = requests.get(url)

        return response
    return duckdb, get_pubchem_id, json, mo, requests


@app.cell
def _(mo):
    data_path = str(mo.notebook_location() / "data" / "jump_pubchem.parquet")
    mapper = mo.sql(
        f"""
        SELECT *,(((from_hex(fingerprint)::BITSTRING)::VARCHAR)[33:-8]) AS bs FROM read_parquet('{data_path}');
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
def _():
    return


@app.cell
def _(duckdb, field, get_pubchem_id, mapper, mo, res, sb):
    if sb.value:
        sim_col = "similarity"
        result = get_pubchem_id(sb.value, field.value)
        fmt = result.json()
        if "Fault" in fmt:
            mo.md("ID not found")
        else:
        #assert "fault" not in fmt, "No CID Found"
            options = [{x["urn"]["label"]:list(x["value"].values())[0] for x in hit["props"] if (x["urn"]["label"]=="InChIKey") or (x["urn"]["label"] == "IUPAC Name" and x["urn"]["name"]=="Preferred") or (x["urn"]["label"]=="Fingerprint") } for hit in fmt["PC_Compounds"]]
            # Convert to bitstring
            fp_hex = options[0]["Fingerprint"]
            bitstring = (bin(int(fp_hex, 16))).zfill(len(fp_hex)*4)[32:-7]

            # Query against JUMP->PubChem mapper
            res = duckdb.sql(
                f"SELECT Metadata_JCP2022,InChIKey,BIT_COUNT(bs::BITSTRING & '{bitstring}') / BIT_COUNT(bs::BITSTRING | '{bitstring}') AS {sim_col} FROM mapper ORDER BY {sim_col} DESC LIMIT 3"
            )
            mo.sql(
            f"""
            SELECT * FROM res
            """
        )
    return bitstring, fmt, fp_hex, options, res, result, sim_col


if __name__ == "__main__":
    app.run()
