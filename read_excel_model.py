import pandas as pd
import re

def read_excel_model():
    species = []
    file = 'CBB_model.xlsx'
    df = pd.read_excel(file, header=None)
    reactions = df.iloc[:, 0].tolist()
    nrxn = len(reactions)
    Si = []  # The row number that is being sparsed
    Sj = []  # The corresponding column number that is being sparsed
    Ss = []  # The matrix S
    rate_inds = [[] for _ in range(nrxn)]  # creates rate_inds list, a 59x1 empty list
    for irxn in range(nrxn):
        rstring = reactions[irxn]
        sep_ind = rstring.find('->')
        lhs = rstring[:sep_ind].strip()
        reactants = [r.strip() for r in lhs.split('+')]
        pattern = r'(?P<stoich>\(?[0-9.]+\)?\s+|)(?P<species>\S+)'
        rate_inds[irxn] = []
        for ireactant, reactant in enumerate(reactants):
            if not reactant:  # Skip empty reactants
                continue
            match = re.match(pattern, reactant)
            if match is None:
                raise ValueError(f"Failed to match reactant '{reactant}' in reaction '{rstring}'")
            stoich_str = match.group('stoich').strip('()').strip()
            species_str = match.group('species')
            stoich = float(stoich_str) if stoich_str else 1.0
            if species_str not in species:
                species.append(species_str)
            ireactant_ind = species.index(species_str)
            Si.append(ireactant_ind)
            Sj.append(irxn)
            Ss.append(-stoich)
            rate_inds[irxn].append(ireactant_ind)
        rhs = rstring[sep_ind + 2:].strip()
        products = [p.strip() for p in rhs.split('+')]
        for iproduct, product in enumerate(products):
            if not product:  # Skip empty products
                continue
            match = re.match(pattern, product)
            if match is None:
                raise ValueError(f"Failed to match product '{product}' in reaction '{rstring}'")
            stoich_str = match.group('stoich').strip('()').strip()
            species_str = match.group('species')
            stoich = float(stoich_str) if stoich_str else 1.0
            if species_str not in species:
                species.append(species_str)
            iproduct_ind = species.index(species_str)
            Si.append(iproduct_ind)
            Sj.append(irxn)
            Ss.append(stoich)
    from scipy.sparse import csr_matrix
    S = csr_matrix((Ss, (Si, Sj)), shape=(len(species), nrxn))
    return species, S.toarray(), rate_inds

if __name__ == "__main__":
    species, S, rate_inds = read_excel_model()
    print("Species:", species)
    print("Sparse Matrix:\n", S)
    print("Rate Indices:", rate_inds)


