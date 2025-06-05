import pdfplumber
import pandas as pd
from cobra.io import read_sbml_model
from bioservices import BiGG

MODEL_PATH = 'iJO1366.xml'
PDF_PATH = '00A-042-Rev-D-Phenotype-MicroArrays-1-10-Plate-Maps.pdf'


def extract_pm1_substrates(pdf_path: str) -> list[str]:
    """Extract substrate names from the PM1 plate map PDF."""
    with pdfplumber.open(pdf_path) as pdf:
        table = pdf.pages[0].extract_table()
    substrates = []
    for row in table:
        for cell in row:
            if not cell:
                continue
            cell = cell.strip()
            if '\n' in cell:
                _well, name = cell.split('\n', 1)
            else:
                name = cell
            name = name.replace('\n', ' ').strip()
            substrates.append(name)
    if substrates and substrates[0].lower().startswith('negative'):
        substrates = substrates[1:]
    return substrates


def normalize_biolog_name(name: str) -> str:
    s = name.lower()
    s = s.replace('α-', '').replace('alpha-', '').replace('β-', '').replace('beta-', '')
    s = s.replace(' ', '_').replace('-', '_').replace(',', '')
    return s


def main():
    model = read_sbml_model(MODEL_PATH)
    pm1_list = extract_pm1_substrates(PDF_PATH)
    bigg = BiGG()
    records = []

    for chem in pm1_list:
        norm_name = normalize_biolog_name(chem)
        try:
            results = bigg.search(norm_name, "metabolites")
        except Exception as e:
            records.append([chem, norm_name, None, None, f'BiGG search failed: {e}'])
            continue
        candidates = [r for r in results if r.get('organism', '').lower() in ('escherichia_coli', '/')]
        hit = None
        if len(candidates) == 1:
            hit = candidates[0]
        elif len(candidates) > 1:
            for r in candidates:
                if r['bigg_id'].lower() == norm_name:
                    hit = r
                    break
            if hit is None:
                hit = candidates[0]
        if not hit:
            records.append([chem, norm_name, None, None, 'not found in BiGG'])
            continue
        base_id = hit['bigg_id']
        if base_id.endswith('_e') or base_id.endswith('_c') or base_id.endswith('_p'):
            met_id = base_id
        else:
            met_id = base_id + '_e'
        try:
            met = model.metabolites.get_by_id(met_id)
        except KeyError:
            records.append([chem, norm_name, base_id, met_id, 'metabolite missing'])
            continue
        ex_reactions = [rxn.id for rxn in met.reactions if rxn.id.startswith('EX_')]
        if ex_reactions:
            records.append([chem, norm_name, base_id, met_id, ex_reactions[0]])
        else:
            records.append([chem, norm_name, base_id, met_id, 'exchange missing'])

    df = pd.DataFrame(records, columns=['Biolog name', 'normalized', 'BiGG base_id', 'Metabolite ID', 'Exchange or status'])
    print(df)
    df.to_csv('pm1_mapping_results.csv', index=False, encoding='utf-8-sig')


if __name__ == '__main__':
    main()
