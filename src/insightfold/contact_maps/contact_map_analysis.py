import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")


@app.cell
def _():
    import io
    import os
    import gemmi
    import marimo as mo
    import requests
    import numpy as np
    import tempfile
    import gzip
    from collections import defaultdict
    import pandas as pd
    from functools import lru_cache
    from concurrent.futures import ThreadPoolExecutor

    return ThreadPoolExecutor, defaultdict, gemmi, lru_cache, mo, requests


@app.cell
def _(mo):
    form = mo.ui.text_area(placeholder="P00519").form()
    mo.vstack([mo.md('# Input Uniprot accession'), form])
    return (form,)


@app.cell
def _(form, get_uniprot_data):
    if form.value is not None:
        uniprot_data = get_uniprot_data(form.value)
        entry_data = get_entry_details(uniprot_data[form.value]['data'])
    return entry_data, uniprot_data


@app.cell
def _(entry_data, form, mo):
    if form.value is not None:
        mo.vstack([mo.md(f'# PDB entries mapping to the protein {form.value}'), 
        mo.ui.table(data=entry_data, pagination=True)])
    return


@app.cell
def _(form, uniprot_data):
    if form.value is not None:
        mutation_data = get_mutation_info(uniprot_data[form.value]['data'])
    return (mutation_data,)


@app.cell
def _(mo):
    get_selected_variant, set_selected_variant = mo.state(None)
    return get_selected_variant, set_selected_variant


@app.cell
def _(form, mo, mutation_data, set_selected_variant):
    def make_view_button(entry_id, uniprot_id, unp_seq_id):
        def on_click(click_count):
            next_count = (click_count or 0) + 1
            set_selected_variant((entry_id, uniprot_id, unp_seq_id, next_count))
            return next_count

        return mo.ui.button(
            value=0,
            on_click=on_click,
            label="View 3D",
            kind="warn",
        )

    mutation_rows = []
    for item in mutation_data:
        row = item.copy()
        row["view"] = make_view_button(
            row["entry"],
            form.value,
            row["uniprot_seq_id"],
        )
        mutation_rows.append(row)
    return (mutation_rows,)


@app.cell
def _(form, get_selected_variant, mo, mutation_rows, visualise_variant):
    selected_variant = get_selected_variant()

    if (
        selected_variant is None
        or form.value is None
        or selected_variant[1] != form.value
    ):
        viewer = mo.md("Click **View 3D** in the table to render the structure viewer.")
    else:
        viewer = visualise_variant(*selected_variant[:3])

    mo.vstack([
        mo.md('# PDB entries with variants'),
        mo.ui.table(data=mutation_rows, pagination=True),
        viewer
    ])

    return


@app.cell
def _(mo):
    mo.md("""
    # Compare variant to wild type
    """)
    return


@app.cell
def _(entry_data, form, mutation_data):
    if form.value is not None:
        variant_entries = sorted({item['entry'] for item in mutation_data})
        wild_entries = sorted({item['entry'] for item in entry_data if item['entry'] not in variant_entries})
    return variant_entries, wild_entries


@app.cell
def _(mo, variant_entries, wild_entries):
    variant_entry_dropdown = mo.ui.dropdown(
        options=variant_entries,
        value=variant_entries[0] if variant_entries else None,
        label="Select a variant entry",
    )
    wild_entry_dropdown = mo.ui.dropdown(
        options=wild_entries,
        value=wild_entries[0] if wild_entries else None,
        label="Select a wild entry",
    )

    mo.hstack([
        variant_entry_dropdown,
        wild_entry_dropdown,
    ])

    return variant_entry_dropdown, wild_entry_dropdown


@app.cell(hide_code=True)
def comapre_wt_to_variant(form, mo, superpose, viz_superpose):
    def comapre_wt_to_variant(wild_entry_id, variant_entry_id, mutation_data):
        if not wild_entry_id or not variant_entry_id:
            return mo.md("Select both a wild-type and variant entry.")
        if not form.value:
            return mo.md("Enter a UniProt accession first.")

        sup_result = superpose(wild_entry_id, variant_entry_id, form.value)
        return viz_superpose(
            wild_entry_id,
            variant_entry_id,
            mutation_data,
            sup_result.transform.vec,
            sup_result.transform.mat,
        )


    return (comapre_wt_to_variant,)


@app.cell
def _(
    comapre_wt_to_variant,
    form,
    mo,
    mutation_data,
    variant_entries,
    variant_entry_dropdown,
    wild_entries,
    wild_entry_dropdown,
):
    def show_superpose():
        if not variant_entries or not wild_entries:
            return mo.md("Need at least one variant entry and one wild-type entry to compare.")
        elif form.value is None:
            return mo.md("Enter a UniProt accession first.")
        elif not wild_entry_dropdown.value or not variant_entry_dropdown.value:
            return mo.md("Select both a wild-type and variant entry to render the aligned view.")
        else:
            return comapre_wt_to_variant(wild_entry_dropdown.value, variant_entry_dropdown.value, mutation_data)
    

    return (show_superpose,)


@app.cell
def _(show_superpose):
    show_superpose()
    return


@app.cell
def _(gemmi, requests):
    def get_entry_assembly(entry_id, assembly_id):
        url = f"https://www.ebi.ac.uk/pdbe/model-server/v1/{entry_id}/assembly?name={assembly_id}&encoding=cif&copy_all_categories=false&download=false"
        response = requests.get(url)
        if response.status_code != 200:
            return
    
        doc = gemmi.cif.read_string(response.content)
        st = gemmi.make_structure_from_block(doc.sole_block())
        return st
                                         

    return (get_entry_assembly,)


@app.cell
def _(get_entry_to_uniprot_map):
    def get_unp_map_polymer(entry_id, model, uniprot_id):
        entry_unp_map = get_entry_to_uniprot_map(entry_id, uniprot_id)
        for item in entry_unp_map:
            if item['auth_asym_id'] in list(map(lambda x: x.name, model)):
                return model[item['auth_asym_id']].get_polymer()

    return (get_unp_map_polymer,)


@app.cell
def _(
    gemmi,
    get_entry_assembly,
    get_preferred_assembly_ids,
    get_unp_map_polymer,
):
    def superpose(entry_1, entry_2, uniprot_id):
        preferred_assemblies = get_preferred_assembly_ids([entry_1, entry_2])
        st_1 = get_entry_assembly(entry_1, str(preferred_assemblies[entry_1]))
        st_2 = get_entry_assembly(entry_2, str(preferred_assemblies[entry_2]))

        polymer_1 = get_unp_map_polymer(entry_1, st_1[0], uniprot_id)
        polymer_2 = get_unp_map_polymer(entry_2, st_2[0], uniprot_id)
        ptype = polymer_1.check_polymer_type()

        sup = gemmi.calculate_superposition(polymer_1, polymer_2, ptype, gemmi.SupSelect.CaP)
    
        return sup

    return (superpose,)


@app.cell
def _(get_entry_to_uniprot_map, is_ligand_binding_site):
    def visualise_variant(entry_id, uniprot_id, unp_seq_id):
        entry_unp_map = get_entry_to_uniprot_map(entry_id, uniprot_id)
        entry_seq_ids = get_uniprot_seq_id_to_entry_seq_id(entry_unp_map, unp_seq_id)
        (binding_ligands, assembly_id) = is_ligand_binding_site(entry_id, uniprot_id, unp_seq_id)
        ligands_to_show = []
        if binding_ligands:
            for lig in binding_ligands:
                ligands_to_show.append({
                    "auth_asym_id": lig[2],
                    "auth_comp_id": lig[1],
                    "auth_seq_id": lig[3],
                    "pdbx_PDB_ins_code": lig[4] or ""
            })

        return viz_variant_with_ligands(entry_id, str(assembly_id), ligands_to_show, entry_seq_ids)

    return (visualise_variant,)


@app.function
def get_uniprot_seq_id_to_entry_seq_id(entry_uniprot_map, uniprot_seq_id):
    entry_seq_ids = []
    for item in entry_uniprot_map:
        if item['unp_start']<=uniprot_seq_id<=item['unp_end']:
            entry_seq_ids.append({
                'auth_asym_id': item['auth_asym_id'],
                'label_seq_id': item['pdb_label_seq_id_start'] + (uniprot_seq_id - item['unp_start'])
            })

    return entry_seq_ids


@app.cell
def _(requests):
    def get_uniprot_data(uniprot_id: str):
        url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot_id}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data

    return (get_uniprot_data,)


@app.function
def get_entry_details(unipdb_data):
    entry_data = []
    for entry in unipdb_data:
        entry_data.append(
            {
                'entry': entry['name'],
                'resolution': entry['additionalData']['resolution'],
                'experiment': entry['additionalData']['experiment'],
                'title': entry['additionalData']['title'],
                'ligand_count': entry['additionalData']['ligandCount']
            }
        )

    return entry_data


@app.function
def get_mutation_info(unipdb_data):
    mutation_data = []
    for entry in unipdb_data:
        for residue in entry['residues']:
            if residue['startIndex'] == residue['endIndex']:
                if ('mutation' in residue):
                    mutation_data.append({
                        'entry': entry['name'],
                        'change': residue['mutationType'],
                        'uniprot_seq_id' : residue['startIndex'],
                        'from_aa' : residue['startCode'],
                        'to_aa': residue['pdbCode']
                    })
    return mutation_data


@app.cell
def _(defaultdict, requests):
    def get_ligand_data(uniprot_id):
        """Get ligand interaction data for a UniProt ID.
        Returns (lig_site, ligands) where lig_site maps ligand accession to binding site seq_ids,
        and ligands maps ligand accession to list of PDB IDs."""
        if not uniprot_id:
            return {}, {}

        ligand_url = f"https://www.ebi.ac.uk/pdbe/api/v2/uniprot/ligands/{uniprot_id}"
        lig_data = None
        lig_response = requests.get(ligand_url)
        if lig_response.status_code == 200:
            lig_data = lig_response.json()

        lig_site = defaultdict(list)
        ligands = defaultdict(list)

        if lig_data:
            for item in lig_data[uniprot_id]:
                for lig_acc, lig_info in item.items():
                    if lig_info['directly_interacts'] == 'true':
                        ligands[lig_acc].extend(lig_info.get('pdbs', []))
                        for site in lig_info.get('binding_sites', []):
                            lig_site[lig_acc].append(site)

        return lig_site, ligands

    return


@app.cell
def _(requests):
    def get_bound_molecules(entry_id):
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_molecules/{entry_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return []

        data = response.json()
        records = []

        for bm in data.get(entry_id, []):
            bm_id = bm.get("bm_id")
            for ligand in bm.get("composition", {}).get("ligands", []):
                pdb_ins_code = ligand.get("author_insertion_code")
                if pdb_ins_code is not None:
                    pdb_ins_code = pdb_ins_code.strip() or None

                records.append({
                    "entry_id": entry_id,
                    "bm_id": bm_id,
                    "auth_asym_id": ligand.get("chain_id"),
                    "auth_seq_id": ligand.get("author_residue_number"),
                    "pdb_ins_code": pdb_ins_code,
                    "auth_comp_id": ligand.get("chem_comp_id"),
                    "entity_id": ligand.get("entity"),
                })

        return records

    return (get_bound_molecules,)


@app.function
def get_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


@app.cell
def _(requests):
    def get_preferred_assembly(entry_id):
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/{entry_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return {}

        data = response.json()
        for assembly in data[entry_id][0]['assemblies']:
            if assembly['preferred']:
                return int(assembly['assembly_id'])




    return (get_preferred_assembly,)


@app.cell
def _(requests):
    def get_preferred_assembly_ids_batch(entry_ids: list[str]):
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/{','.join(entry_ids)}"
        response = requests.get(url)
        if response.status_code != 200:
            return {}

        data = response.json()
        pref_assembly_ids = {}
        for entry_id in data:
            for assembly in data[entry_id][0]['assemblies']:
                if assembly['preferred']:
                    pref_assembly_ids[entry_id] = int(assembly['assembly_id'])
                    break

        return pref_assembly_ids

    return (get_preferred_assembly_ids_batch,)


@app.cell
def _(ThreadPoolExecutor, get_preferred_assembly_ids_batch):
    def get_preferred_assembly_ids(entry_ids: list[str]):
        """Get preferred assembly IDs for a list of PDB entries using batched API calls."""
        chunk_size = 50
        chunks = list(get_chunks(entry_ids, chunk_size))
        results = {}
        with ThreadPoolExecutor() as exe:
            for result in exe.map(get_preferred_assembly_ids_batch, chunks):
                results.update(result)
        return results

    return (get_preferred_assembly_ids,)


@app.cell
def _(requests):
    def get_ligand_interactions(pdb_id, auth_asym_id, auth_seq_id):
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_ligand_interactions/{pdb_id}/{auth_asym_id}/{auth_seq_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return []

        data = response.json()
        records = []

        for pdb_entry, entry_data in data.items():
            for ligand_info in entry_data:
                ligand = ligand_info.get("ligand", {})
                interactions = ligand_info.get("interactions", [])

                for interaction in interactions:
                    ligand_atoms = interaction.get("ligand_atoms", [])
                    interaction_type = interaction.get("interaction_type", "")
                    interaction_details = interaction.get("interaction_details", [])
                    distance = interaction.get("distance")
                    end = interaction.get("end", {})
                    target_atoms = end.get("atom_names", [])

                    base = {
                        "pdb_id": pdb_entry,
                        "lig_auth_asym_id": ligand.get("chain_id"),
                        "lig_auth_seq_id": ligand.get("author_residue_number"),
                        "lig_comp_id": ligand.get("chem_comp_id"),
                        "interaction_type": interaction_type,
                        "distance": distance,
                        "target_auth_asym_id": end.get("chain_id"),
                        "target_auth_seq_id": end.get("author_residue_number"),
                        "target_comp_id": end.get("chem_comp_id"),
                    }

                    for lig_atom in ligand_atoms:
                        for target_atom in target_atoms:
                            for detail in interaction_details:
                                records.append({
                                    **base,
                                    "lig_atom": lig_atom,
                                    "target_atom": target_atom,
                                    "interaction_type": detail,
                                })

        return records

    return


@app.cell
def _(lru_cache, requests):
    @lru_cache()
    def get_entry_to_uniprot_map(entry_id, uniprot_id):
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return []

        data = response.json()
        mappings = data.get(entry_id, {}).get("UniProt", {}).get(uniprot_id, {})
        records = []

        for mapping in mappings.get("mappings", []):
            records.append({
                "entity_id": mapping.get("entity_id"),
                "auth_asym_id": mapping.get("chain_id"),
                "unp_start": mapping.get("unp_start"),
                "unp_end": mapping.get("unp_end"),
                "pdb_auth_seq_id_start": mapping.get("start", {}).get("author_residue_number"),
                "pdb_auth_seq_id_end": mapping.get("end", {}).get("author_residue_number"),
                "pdb_label_seq_id_start": mapping.get("start", {}).get("residue_number"),
                "pdb_label_seq_id_end": mapping.get("end", {}).get("residue_number"),
                "identity": mapping.get("identity"),
                "coverage": mapping.get("coverage"),
            })

        return records

    return (get_entry_to_uniprot_map,)


@app.cell
def _(gemmi, requests):
    def get_interacting_residues(entry_id, assembly_id, entity_id, auth_asym_id, auth_comp_id, auth_seq_id, pdb_ins_code=None):
        url = f"https://www.ebi.ac.uk/pdbe/model-server/v1/{entry_id}/residueInteraction?label_entity_id={entity_id}&auth_asym_id={auth_asym_id}&auth_comp_id={auth_comp_id}&auth_seq_id={auth_seq_id}&radius=5&assembly_name={assembly_id}&encoding=cif&copy_all_categories=false&download=false"
        response = requests.get(url)
        if response.status_code != 200:
            return []

        intx_residues = []
        block = gemmi.cif.read_string(response.content).sole_block()
        table = block.find(['_atom_site.group_PDB', '_atom_site.label_asym_id', '_atom_site.label_seq_id', '_atom_site.label_comp_id', '_atom_site.pdbx_sifts_xref_db_name', '_atom_site.pdbx_sifts_xref_db_acc', '_atom_site.pdbx_sifts_xref_db_num'])

        current_seq_id = None
        for row in table:
            if row[0] == 'ATOM':
                if current_seq_id != f"{row[2]}_{row[3]}":
                    intx_residues.append(dict(zip(["label_asym_id","label_seq_id","label_comp_id","sifts_db", "sifts_db_acc", "sifts_residue_num"], [row[i] for i in range(1, table.width())])))
                    current_seq_id = f"{row[2]}_{row[3]}"

        return intx_residues

    return (get_interacting_residues,)


@app.cell
def _(
    get_bound_molecules,
    get_interacting_residues,
    get_preferred_assembly,
    lru_cache,
):
    @lru_cache()
    def get_entry_binding_sites(entry_id: str):
        pref_assembly_id = get_preferred_assembly(entry_id)
        chem_comps = get_bound_molecules(entry_id)

        entry_intx_residues = {}
        for chem_comp in chem_comps:
            entry_id = chem_comp["entry_id"]
            if pref_assembly_id is None:
                continue

            intx_residues = get_interacting_residues(
                entry_id,
                pref_assembly_id,
                chem_comp["entity_id"],
                chem_comp["auth_asym_id"],
                chem_comp["auth_comp_id"],
                chem_comp["auth_seq_id"],
                chem_comp["pdb_ins_code"],
            )

            ligand_instance = (
                chem_comp["entity_id"],
                chem_comp["auth_comp_id"],
                chem_comp["auth_asym_id"],
                chem_comp["auth_seq_id"],
                chem_comp["pdb_ins_code"],
            )
            entry_intx_residues.setdefault(entry_id, {})[ligand_instance] = intx_residues

        return (entry_intx_residues, pref_assembly_id)

    return (get_entry_binding_sites,)


@app.cell
def _(
    ThreadPoolExecutor,
    get_bound_molecules,
    get_interacting_residues,
    get_preferred_assembly_ids,
):
    def get_binding_sites(entry_ids: list[str]):
        pref_assembly_ids = get_preferred_assembly_ids(entry_ids)
        chem_comps = []
        with ThreadPoolExecutor() as executor:
            for records in executor.map(get_bound_molecules, entry_ids):
                chem_comps.extend(records)

        entry_intx_residues = {}
        for chem_comp in chem_comps:
            entry_id = chem_comp["entry_id"]
            assembly_id = pref_assembly_ids.get(entry_id)
            if assembly_id is None:
                continue

            intx_residues = get_interacting_residues(
                entry_id,
                assembly_id,
                chem_comp["entity_id"],
                chem_comp["auth_asym_id"],
                chem_comp["auth_comp_id"],
                chem_comp["auth_seq_id"],
                chem_comp["pdb_ins_code"],
            )

            ligand_instance = (
                chem_comp["entity_id"],
                chem_comp["auth_comp_id"],
                chem_comp["auth_asym_id"],
                chem_comp["auth_seq_id"],
                chem_comp["pdb_ins_code"],
            )
            entry_intx_residues.setdefault(entry_id, {})[ligand_instance] = intx_residues

        return entry_intx_residues

    return


@app.cell
def _(defaultdict, get_entry_binding_sites):
    def is_ligand_binding_site(entry_id, uniprot_id, uniprot_seq_id):
        mutated_bs_sites = defaultdict(list)
        (bs, assembly_id) = get_entry_binding_sites(entry_id)
        if not bs:
            return (None, None)
        for ligand, residues in bs[entry_id].items():
            for residue in residues:
                if residue['sifts_db_acc'] == uniprot_id:
                    if int(residue['sifts_residue_num']) == uniprot_seq_id:
                        mutated_bs_sites[ligand].append(residue)


        return (mutated_bs_sites, assembly_id)

    return (is_ligand_binding_site,)


@app.function
def viz_variant_with_ligands(entry_id, assembly_id, ligands, residues, write_output=False):
    """Build a Mol* view with selected ligands and residues in ball-and-stick."""
    mo = __import__("marimo")
    mvs = __import__("molviewspec")
    mvs_builder_module = __import__("molviewspec.builder", fromlist=["create_builder"])

    ligand_required_keys = {"auth_asym_id", "auth_comp_id", "auth_seq_id", "pdbx_PDB_ins_code"}
    residue_required_keys = {"auth_asym_id", "label_seq_id"}
    pdb_id = str(entry_id).strip().lower()
    if not pdb_id:
        raise ValueError("entry_id must be a non-empty PDB entry id")
    def normalise_ligand_selector(ligand):
        missing = ligand_required_keys - set(ligand)
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValueError(f"ligands entry is missing: {missing_text}")
        return {
            "auth_asym_id": ligand["auth_asym_id"],
            "auth_comp_id": ligand["auth_comp_id"],
            "auth_seq_id": int(ligand["auth_seq_id"]),
            "pdbx_PDB_ins_code": ligand["pdbx_PDB_ins_code"],
        }

    def normalise_residue_selector(residue):
        missing = residue_required_keys - set(residue)
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValueError(f"residues entry is missing: {missing_text}")
        return {
            "auth_asym_id": residue["auth_asym_id"],
            "label_seq_id": int(residue["label_seq_id"]),
        }

    ligand_selectors = [normalise_ligand_selector(ligand) for ligand in ligands]
    residue_selectors = [normalise_residue_selector(residue) for residue in residues]

    builder = mvs_builder_module.create_builder()
    structure = (
        builder
        .download(url=f"https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id}.bcif")
        .parse(format="bcif")
        .assembly_structure(assembly_id=assembly_id)
    )

    polymer = structure.component(selector="polymer")
    polymer.representation(type="cartoon").color(color="#6f8fa6")

    def add_ball_and_stick(selector, color, label, tooltip):
        component = structure.component(selector=selector)
        component.representation(type="ball_and_stick", size_factor=1).color(color=color)
        component.label(text=label)
        component.tooltip(text=tooltip)
        return component

    for selector in ligand_selectors:
        insertion_code = selector["pdbx_PDB_ins_code"]
        label = f"{selector['auth_comp_id']} {selector['auth_asym_id']} {selector['auth_seq_id']}"
        if insertion_code not in (None, " ", "?", "."):
            label = f"{label}{insertion_code}"
        add_ball_and_stick(
            selector=selector,
            color="#1f77b4",
            label=label,
            tooltip=(
                f"Ligand: {selector['auth_comp_id']}, "
                f"auth_asym_id {selector['auth_asym_id']}, "
                f"auth_seq_id {selector['auth_seq_id']}, "
                f"pdbx_PDB_ins_code {selector['pdbx_PDB_ins_code']}"
            ),
        )

    for selector in residue_selectors:
        label = f"{selector['auth_asym_id']} {selector['label_seq_id']}"
        add_ball_and_stick(
            selector=selector,
            color="#d62728",
            label=label,
            tooltip=(
                f"auth_asym_id {selector['auth_asym_id']}, "
                f"label_seq_id {selector['label_seq_id']}"
            ),
        )

    focus_selectors = ligand_selectors + residue_selectors
    if focus_selectors:
        structure.component(selector=focus_selectors).focus(radius_factor=1)

    builder.canvas(background_color="white")
    state = builder.get_state(title=f"{pdb_id.upper()} selected ligands and residues")
    mvsj = mvs.MVSJ(data=state)

    if write_output:
        output_filename = f"{pdb_id}_selected_ligands_residues.mvsj"
        mvsj.dump(output_filename)

    html = mvs.molstar_html(mvsj, molstar_version="latest")
    return mo.iframe(html, height="650px")


@app.cell(hide_code=True)
def viz_superpose(
    form,
    get_entry_to_uniprot_map,
    get_preferred_assembly_ids,
    is_ligand_binding_site,
):
    def viz_superpose(entry_1, entry_2, mutation_data, transformation_matrix, rotation_matrix, write_output=False):
        """Build a Mol* view of two superposed structures using MolViewSpec.

        The function accepts either of these call styles:
        - ``viz_superpose(entry_1, entry_2, mutation_data, sup.transform.vec, sup.transform.mat)``
        - ``viz_superpose(entry_1, entry_2, mutation_data, sup.transform.mat, sup.transform.vec)``

        The first entry is treated as wild type. Residues corresponding to the
        variant UniProt positions are highlighted on ``entry_1`` via UniProt
        mapping, while variant residues and nearby ligands are highlighted on
        ``entry_2``.
        """
        mo = __import__("marimo")
        mvs = __import__("molviewspec")
        np = __import__("numpy")
        mvs_builder_module = __import__("molviewspec.builder", fromlist=["create_builder"])

        pdb_id_1 = str(entry_1).strip().lower()
        pdb_id_2 = str(entry_2).strip().lower()
        if not pdb_id_1 or not pdb_id_2:
            raise ValueError("entry ids must be non-empty PDB identifiers")

        preferred_assemblies = get_preferred_assembly_ids([pdb_id_1, pdb_id_2])
        if pdb_id_1 not in preferred_assemblies or pdb_id_2 not in preferred_assemblies:
            raise ValueError("Could not resolve preferred assembly ids for both entries")

        uniprot_id = getattr(form, "value", None)
        if not uniprot_id:
            raise ValueError("A UniProt accession is required in form.value to locate variant residues")

        polymer_1_color = "#4c78a8"
        polymer_2_color = "#e4572e"

        def _to_plain_array(value):
            if hasattr(value, "tolist"):
                value = value.tolist()
            return np.asarray(value, dtype=float)

        def _normalise_rotation(value):
            rotation = _to_plain_array(value)
            if rotation.shape == (3, 3):
                return rotation.flatten(order="F").tolist()
            if rotation.shape == (9,):
                return rotation.tolist()
            raise ValueError("rotation matrix must be a 3x3 matrix or a flat length-9 vector")

        def _normalise_translation(value):
            translation = _to_plain_array(value)
            if translation.shape == (3,):
                return translation.tolist()
            if translation.shape == (1, 3):
                return translation.reshape(3).tolist()
            raise ValueError("translation vector must be a length-3 vector")

        def _dedupe_dicts(items):
            seen = set()
            output = []
            for item in items:
                key = tuple(sorted(item.items()))
                if key in seen:
                    continue
                seen.add(key)
                output.append(item)
            return output

        def _build_annotation_map(entry_unp_map, rows, residue_role):
            selectors = []
            annotations = {}
            for row in rows:
                seq_id = int(row["uniprot_seq_id"])
                mapped_selectors = get_uniprot_seq_id_to_entry_seq_id(entry_unp_map, seq_id)
                for selector in mapped_selectors:
                    norm_selector = {
                        "auth_asym_id": selector["auth_asym_id"],
                        "label_seq_id": int(selector["label_seq_id"]),
                    }
                    selectors.append(norm_selector)
                    key = (norm_selector["auth_asym_id"], norm_selector["label_seq_id"])
                    if residue_role == "wild_type":
                        annotations[key] = {
                            "label": f"WT {row['from_aa']}{seq_id}",
                            "tooltip": (
                                f"Wild-type residue corresponding to variant site in {pdb_id_1.upper()}: "
                                f"{row['from_aa']} {seq_id}"
                            ),
                        }
                    else:
                        annotations[key] = {
                            "label": f"{row['from_aa']}{seq_id}{row['to_aa']}",
                            "tooltip": (
                                f"Variant in {pdb_id_2.upper()}: {row['from_aa']} {seq_id} -> {row['to_aa']} "
                                f"({row['change']})"
                            ),
                        }
            return _dedupe_dicts(selectors), annotations

        first = _to_plain_array(transformation_matrix)
        second = None if rotation_matrix is None else _to_plain_array(rotation_matrix)

        if first.shape == (4, 4):
            transform_kwargs = {"matrix": first.flatten(order="F").tolist()}
        elif first.shape == (16,):
            transform_kwargs = {"matrix": first.tolist()}
        else:
            rotation_shapes = {(3, 3), (9,)}
            translation_shapes = {(3,), (1, 3)}

            if first.shape in translation_shapes and second is not None and second.shape in rotation_shapes:
                translation_value = transformation_matrix
                rotation_value = rotation_matrix
            elif first.shape in rotation_shapes and second is not None and second.shape in translation_shapes:
                translation_value = rotation_matrix
                rotation_value = transformation_matrix
            else:
                raise ValueError(
                    "Expected a rotation matrix and translation vector from superpose(); "
                    "accepted orders are (vec, mat) or (mat, vec)"
                )

            transform_kwargs = {
                "translation": _normalise_translation(translation_value),
                "rotation": _normalise_rotation(rotation_value),
            }

        variant_rows = [row for row in mutation_data if str(row.get("entry", "")).strip().lower() == pdb_id_2]

        entry_1_unp_map = get_entry_to_uniprot_map(pdb_id_1, uniprot_id)
        entry_2_unp_map = get_entry_to_uniprot_map(pdb_id_2, uniprot_id)

        wt_residue_selectors, wt_annotations = _build_annotation_map(entry_1_unp_map, variant_rows, "wild_type")
        variant_residue_selectors, variant_annotations = _build_annotation_map(entry_2_unp_map, variant_rows, "variant")

        ligand_selectors = []
        for row in variant_rows:
            seq_id = int(row["uniprot_seq_id"])
            binding_ligands, _ = is_ligand_binding_site(pdb_id_2, uniprot_id, seq_id)
            if binding_ligands:
                for ligand in binding_ligands:
                    ligand_selectors.append(
                        {
                            "auth_asym_id": ligand[2],
                            "auth_comp_id": ligand[1],
                            "auth_seq_id": int(ligand[3]),
                            "pdbx_PDB_ins_code": ligand[4] or "",
                        }
                    )
        ligand_selectors = _dedupe_dicts(ligand_selectors)

        builder = mvs_builder_module.create_builder()

        structure_1 = (
            builder
            .download(url=f"https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id_1}.bcif")
            .parse(format="bcif")
            .assembly_structure(assembly_id=str(preferred_assemblies[pdb_id_1]))
        )
        structure_2 = (
            builder
            .download(url=f"https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id_2}.bcif")
            .parse(format="bcif")
            .assembly_structure(assembly_id=str(preferred_assemblies[pdb_id_2]))
            .transform(**transform_kwargs)
        )

        polymer_1 = structure_1.component(selector="polymer")
        polymer_1.representation(type="cartoon").color(color=polymer_1_color)

        polymer_2 = structure_2.component(selector="polymer")
        polymer_2.representation(type="cartoon").color(color=polymer_2_color)

        for selector in wt_residue_selectors:
            key = (selector["auth_asym_id"], selector["label_seq_id"])
            annotation = wt_annotations[key]
            component = structure_1.component(selector=selector)
            component.representation(type="ball_and_stick", size_factor=1).color(color=polymer_1_color)
            component.tooltip(text=annotation["tooltip"])

        for selector in ligand_selectors:
            insertion_code = selector["pdbx_PDB_ins_code"]
            label = f"{selector['auth_comp_id']} {selector['auth_asym_id']} {selector['auth_seq_id']}"
            if insertion_code not in (None, "", " ", "?", "."):
                label = f"{label}{insertion_code}"

            component = structure_2.component(selector=selector)
            component.representation(type="ball_and_stick", size_factor=1).color(color=polymer_2_color)
            component.label(text=label)
            component.tooltip(
                text=(
                    f"Ligand near variant in {pdb_id_2.upper()}: {selector['auth_comp_id']} "
                    f"{selector['auth_asym_id']} {selector['auth_seq_id']}"
                )
            )

        for selector in variant_residue_selectors:
            key = (selector["auth_asym_id"], selector["label_seq_id"])
            annotation = variant_annotations[key]
            component = structure_2.component(selector=selector)
            component.representation(type="ball_and_stick", size_factor=1).color(color=polymer_2_color)
            component.tooltip(text=annotation["tooltip"])

        focus_selectors = wt_residue_selectors + variant_residue_selectors + ligand_selectors
        if focus_selectors:
            structure_2.component(selector=variant_residue_selectors + ligand_selectors).focus(radius_factor=1.3)
        else:
            polymer_1.focus(radius_factor=1.2)

        builder.canvas(background_color="white")
        state = builder.get_state(title=f"{pdb_id_1.upper()} vs {pdb_id_2.upper()} superposition")
        mvsj = mvs.MVSJ(data=state)

        if write_output:
            output_filename = f"{pdb_id_1}_{pdb_id_2}_superposed.mvsj"
            mvsj.dump(output_filename)

        html = mvs.molstar_html(mvsj, molstar_version="latest")
        return mo.iframe(html, height="650px")


    viz_supoerpose = viz_superpose

    return (viz_superpose,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
