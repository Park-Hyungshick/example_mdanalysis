from openmm.app import ForceField
import MDAnalysis as mda
import sys

def QfromFF(universe, forcefield):
    residues = universe.residues.resnames
    resnames=forcefield._templates.keys()
    name2charges=dict()

    # Checking 1
    assert len(list(set(residues))) <= len(resnames)

    residue_name_map = dict()
    
    for res_u in list(set(residues)):
        for res_ff in resnames:
            if (res_u in res_ff) or (res_ff in res_u) :
                print(res_u, res_ff)
                
                if( res_u in residue_name_map ): 
                    print("ERROR: Residue names are crashed!")
                    sys.exit()
                else:
                    residue_name_map[res_u] = res_ff
                    
                break
    print(residue_name_map)
    try:
        inverse_map = { v:k for k,v in residue_name_map.items() }
    except:
        print("Resdue name 1-to-1 correspondancy error!")
        sys.exit()

    # Set atom attribute
    universe.add_TopologyAttr('charges')

    # Collect charge info
    for res in resnames:
        charges=[]
        for atom in forcefield._templates[res].atoms :
            charges.append(atom.parameters['charge'])
        try:
            name2charges[ inverse_map[res] ] = charges
        except:
            print(f"{res} residue charge mapping is skipped!")

    total_charges = []
    for res in residues :
        total_charges += name2charges[res]

    universe.atoms.__setattr__('charges',total_charges)
    del total_charges

    return 
