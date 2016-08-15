#!/usr/bin/python

from argparse import ArgumentParser
from sys import stderr, stdout

def rechain( pdb ):
    ''' return 
    '''
    chainid_idx = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R',
                  'S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j',
                  'k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','1','2',
                  '3','4','5','6','7','8','9']
    xyzDict_    = {}  # { rsn: { atom:xyz }}
    pdblineDict    = {}  # { rsn: { atom:xyz }}
    chainlist = []
    chainmap = {}
    with open( pdb, "r" ) as f:
        for pdbline in f:
            if pdbline.startswith("ATOM"):
                rsn      = int( pdbline[22:26]) # res number
                res      = pdbline[17:20].strip() # res name 
                atom     = pdbline[12:16].strip() # atom name
                chain_id = pdbline[21] # chain id

                # when transition trigger a new chainmap
                if not chainlist:
                    #print "initiate", chain_id
                    chainmap[ chain_id ] = chain_id
                    chainid_idx.remove( chain_id )
                    chainlist.append( chain_id )
                    prev_chain_id = chain_id

                if chain_id != prev_chain_id:
                    #print "trans", prev_chain_id, chain_id
                    if chain_id in chainlist:
                        #print chain_id, "in", chainlist,
                        x = chainid_idx[0] # get a new chainid
                        #print ", map", chain_id, "to", x
                        chainid_idx.remove( x )
                        chainmap[ chain_id ] = x
                        chainlist.append(x)
                    else:
                        chainmap[ chain_id ] = chain_id
                        chainid_idx.remove( chain_id )
                        chainlist.append(chain_id)

                #xyz = map(float, [ pdbline[30:38], pdbline[38:46], pdbline[46:54] ])

                new_chain_id = chainmap[ chain_id ]
                new_pdbline = pdbline[:21] + "%s" % new_chain_id + pdbline[22:] 
                #print chain_id, new_chain_id
                if new_chain_id not in pdblineDict.keys():
                    pdblineDict[ new_chain_id ] = { rsn: new_pdbline } 
                    #print "new_chain_id", new_chain_id, xyzDict_[ new_chain_id ].keys()
                elif rsn not in pdblineDict[ new_chain_id ].keys():
                    #print xyzDict_[ new_chain_id ]
                    pdblineDict[ new_chain_id ][ rsn ] = new_pdbline
                else:
                    pdblineDict[ new_chain_id ][ rsn ] += new_pdbline

                prev_chain_id = chain_id
            
    return pdblineDict



if __name__=="__main__":
    parser = ArgumentParser()
    #parser.add_argument("-a", "--chain1", nargs="+", required=True,  help="")
    #parser.add_argument("-b", "--chain2", nargs="+", required=True,  help="")
    parser.add_argument("-p", "--pdb", required=True,  help="")
    args = parser.parse_args()
    dict = rechain( args.pdb )
    for chain in dict.keys():
        stderr.write(chain+"\n")
        for rsn in sorted( dict[chain].keys() ):
            stdout.write(dict[ chain ][rsn])
        stdout.write("TER\n")
