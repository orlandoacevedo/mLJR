#!/usr/bin/env python3

import os
import math
import argparse

__version__ = 0.50

# Modified Lydersen-Joback-Reid method
#
# Reference:
#    Paper 1:  Mirza et al. J. Chem. Eng. Data 2015, 60, 1844−1854.
#    Paper 2:  Valderrama et al. Ind. Eng. Chem. Res., 2008, 47, 8416-8422.


# data in an array
# symbol    deltaTbM    deltaTM    deltaPM   deltaVM      nickName(as many as you can)
# 
# Note: for each entry, both symbol and nickNames will be used to identify inputs

DATA_LJR_NO_RING = [
['-CH3'  ,    23.58,    0.0275,    0.3031,    66.81,     'methyl'     ],
['-CH2-' ,    22.88,    0.0159,    0.2165,    57.11,                  ],
['>CH-'  ,    21.74,    0.0002,    0.1140,    45.70,                  ],
['>C<'   ,    18.18,   -0.0206,    0.0539,    21.78,     '[>C-]-'     ],
['=CH2'  ,    24.96,    0.0170,    0.2493,    60.37,                  ],
['=CH-'  ,    18.25,    0.0182,    0.1866,    49.92,                  ],
['=C<'   ,    24.14,   -0.0003,    0.0832,    34.90,                  ],
['=C='   ,    26.15,   -0.0029,    0.0934,    33.85,                  ],
['@CH'   ,     None,    0.0078,    0.1429,    43.97,                  ],   # three bonds
['@C-'   ,     None,    0.0078,    0.1429,    43.97,                  ],   # three bonds
['-OH'   ,    92.88,    0.0723,    0.1343,    30.40,     'alcohol'    ],
['-O-'   ,    22.42,    0.0051,    0.1300,    15.61,  '[-O]-', 'ether'    ],
['>C=O'  ,    94.97,    0.0247,    0.2341,    69.76,     'ketone'     ],
['-CHO'  ,    72.24,    0.0294,    0.3128,    77.46,     'aldehyde'   ],
['-COOH' ,   169.06,    0.0853,    0.4537,    88.60,      'acid'      ],
['-COO-' ,    81.10,    0.0377,    0.4139,    84.76,      'ester'     ],
['HCOO-' ,     None,    0.0360,    0.4752,    97.77,     'formate'    ],
['=O'    ,   -10.50,    0.0273,    0.2042,    44.03,      'oxide'     ],
['-NH2'  ,    73.23,    0.0364,    0.1692,    49.10,      'amine'     ],
['>NH'   ,    50.17,    0.0119,    0.0322,    78.96,                  ],
['>N-'   ,    11.74,   -0.0028,    0.0304,    26.70,   '[>N<]+', 'Ammonium' ],
['-N='   ,    74.60,    0.0172,    0.1541,    45.54,                  ],
['-CN'   ,   125.66,    0.0506,    0.3697,    89.32,    'nitrile'     ],
['-NO2'  ,   152.54,    0.0448,    0.4529,   123.62,    'nitrate'     ],
['-F'    ,    -0.03,    0.0228,    0.2912,    31.47,   '[F-]',  'F-',  'Fluoro'  ],
['-Cl'   ,    38.13,    0.0188,    0.3738,    62.08,   '[Cl-]', 'Cl-', 'Chloro'  ],
['-Br'   ,    66.86,    0.0124,    0.5799,    76.60,   '[Br-]', 'Br-', 'Bromo'   ],
['-I'    ,    93.84,    0.0148,    0.9174,   100.79,   '[I-]',  'I-',  'Iodo'    ],

['-B'    ,   -24.56,    0.0352,    0.0348,    22.45,      'Boron'     ],
['-P'    ,    34.86,   -0.0084,    0.1776,    67.01,    'Phosphoro'   ],
['-SO2'  ,   147.24,   -0.0563,   -0.0606,   112.19,    'sulfoxide'   ],
]


DATA_LJR_RING = [
['-CH2-' ,    27.15,    0.0116,    0.1982,    51.64,                  ],
['>CH-'  ,    21.78,    0.0081,    0.1773,    30.56,                  ],
['=CH-'  ,    26.73,    0.0114,    0.1693,    42.55,                  ],
['>C<'   ,    21.32,   -0.0180,    0.0139,    17.62,                  ],
['=C<'   ,    31.01,    0.0051,    0.0955,    31.28,                  ],
['-O-'   ,    31.22,    0.0138,    0.1371,    17.41,                  ],
['-OH'   ,    76.34,    0.0291,    0.0493,   -17.44,      'phenol'    ],
['>C=O'  ,    94.97,    0.0343,    0.2751,    59.32,                  ],
['>NH'   ,    52.82,    0.0244,    0.0724,    27.61,                  ],
['>N-'   ,    68.16,    0.0063,    0.0538,    25.17,                  ],  #caution, difference, on paper 2
['-N='   ,    57.55,   -0.0011,    0.0559,    42.15,                  ],
]



# a template for Choline Chloride   --  Purity
TEMPLATE_CCL = {
    'type'      : 'purity',
    'name'      : 'Choline Chloride',
    'm'         :  139.62,
    's'         : 'C5H14NO Cl',

    'no-ring' : {
        '-CH3'  : 3,
        '-CH2-' : 2,
        '-OH'   : 1,
        '>N-'   : 1,
        '-Cl'   : 1
    }
}

# a template for Choline Chloride Glycerol   --  mixture
TEMPLATE_CCG = {
    'type'    : 'mixture',
    'name'    : 'Choline Chloride Glycerol',
    'ratio'   : '1:2',
    'm1'      : 139.62,
    's1'      : 'C5H14NO Cl',
    'm2'      : 92.09,
    's2'      : 'C3H8O3',
    't-ref'   :  298.15,
    'q-ref'   :  56.0,
    't'       :  425,

    # Choline Chloride -- molecule 1
    'no-ring-1' : {
        '-CH3'  : 3,
        '-CH2-' : 2,
        '-OH'   : 1,
        '>N-'   : 1,
        '-Cl'   : 1,
    },

    # Glycerol -- molecule 2
    'no-ring-2' : {
        '-CH2-' : 2,
        '-OH'   : 3,
        '>CH-'  : 1,
    }
}



# The following is the predefined info for help message
def func_show_template(template,file=False):
    """This function is used to show help message either to screen or write a template when file is
       set to True. Note: the new generated template file will always NOT overwrite any files"""
    def fout(key):
        txt = ''
        for i,j in key.items(): txt += '    {:7} :  {:>3}\n'.format(i,j)
        return txt

    info = '# Template to calculate Critical Properties either for purity or mixture\n' + \
           '#\n' + \
           '# Free formats. Char < # > is used for comments\n' + \
           '# Case-insensitive, number of spaces or quotes does not matter\n' + \
           '#\n' + \
           '# Note: only keywords < type >; < mark >; < ratio >; < M > < S > & < no-ring > & < ring >\n' + \
           '#       (for purity) OR < M1 > < S1 > & < M2 > < S2 > & < no-ring-1 > & < no-ring-2 > &\n' + \
           '#       < ring-1 > & < ring-2 > (for mixture); < T-ref > < Q-ref > < T > can be used\n' + \
           '#\n' + \
           '# Note: keyword < name > will be ignored\n' + \
           '# Note: keyword < M > is precendent of < S >, they both are used to get molecular weight\n' + \
           '# Note: keyword < S > structrue in a format: (atomType & number) can be separated by "-"\n' + \
           '#       OR blank space. Number of "-" and spaces does not matter\n' + \
           '#            Example 1:   C2H3O4N5\n' + \
           '#            Example 2:   C2 - H3 - O4 - N5\n' + \
           '#            Example 3:   C2   H3   O4   N5\n' + \
           '#            Example 4:   C2   H3 - O4 - N5\n' + \
           '# Note: keyword < t-ref >, < q-ref >, < t > is used to calculate Surface Tension\n' + \
           '#       < t-ref >  :  reference temperature\n' + \
           '#       < q-ref >  :  reference Surface Tension at T-ref\n' + \
           '#         < t >    :  temperature where to calculate Surface Tension\n' + \
           '#\n' + \
           '# All < mark > blocks will be auto assigned to different molecules\n' + \
           '# Inside it, only the defined "Group Name" can be used along with their numbers\n' + \
           '# To check all available "Group Names", please specify "--avail-group-name" or "-g"\n' + \
           '#\n' + \
           '# The number of group can be omitted if it is one\n' + \
           '# For every group, they have to be input line by line\n' + \
           '# Symbol @ means three bonds\n' + \
           '#\n' + \
           '# Specially, for < type >, "mixture" & "mix" & "m" OR "purity" & "pure" & "p" is equal\n\n\n'

    if 'type' in template: info += 'type : ' + template['type'] + '\n\n'
    if 'name' in template: info += 'name : ' + template['name'] + '\n\n'
    if 'ratio' in template: info += 'ratio : ' + template['ratio'] + '\n\n'
    if 'm' in template: info += 'm : ' + str(template['m']) + '\n\n'
    if 'm1' in template: info += 'm1 : ' + str(template['m1']) + '\n\n'
    if 's1' in template: info += 's1 : ' + template['s1'] + '\n\n'
    if 's' in template: info += 's : ' + template['s'] + '\n\n'
    if 'm2' in template: info += 'm2 : ' + str(template['m2']) + '\n\n'
    if 's2' in template: info += 's2 : ' + template['s2'] + '\n\n'
    if 't-ref' in template: info += 't-ref : ' + str(template['t-ref']) + '\n\n'
    if 'q-ref' in template: info += 'q-ref : ' + str(template['q-ref']) + '\n\n'
    if 't' in template: info += 't : ' + str(template['t']) + '\n\n'

    if 'ring' in template: info += 'mark : ring\n\n' + fout(template['ring']) + '\n\n'
    if 'no-ring' in template: info += 'mark : no-ring\n\n' + fout(template['no-ring']) + '\n\n'

    if 'ring-1' in template: info += 'mark : ring-1\n\n' + fout(template['ring-1']) + '\n\n'
    if 'no-ring-1' in template: info += 'mark : no-ring-1\n\n' + fout(template['no-ring-1']) + '\n\n'
    
    if 'ring-2' in template: info += 'mark : ring-2\n\n' + fout(template['ring-2']) + '\n\n'
    if 'no-ring-2' in template: info += 'mark : no-ring-2\n\n' + fout(template['no-ring-2']) + '\n\n'

    if file:
        fname = "template_mLJR.txt"
        cnt = 1
        while True:
            if os.path.isfile(fname):
                fname = 'template_mLJR_{:}.txt'.format(cnt)
            else:
                break
            cnt += 1
        with open(fname,'wt') as f: f.write(info)
        print('Note: template is saved to file < {:} >'.format(fname))
    else:
        print(info)



def func_proDATA(inlist):
    """
    Pre-process input DATA list
    
    Parameter:
        inlist:  2D nested list
        sublist in an array
        [ symbol    deltaTbM    deltaTM    deltaPM   deltaVM      nickName(as many as you can) ]
    
    Return:
        3D nested list

        in a format,

        [
            [ [ allNames ], [ deltaTbM    deltaTM    deltaPM   deltaVM ] ],
            [ [ allNames ], [ deltaTbM    deltaTM    deltaPM   deltaVM ] ],
            ...
        ]
    """
    # avoid any predefined-error
    prolist = []
    bo = False
    for i in inlist:
        if isinstance(i,list) and len(i) >= 5:
            data = []
            for j in i[1:5]:
                if j is None or isinstance(j,(int,float)):
                    data.append(j)
                else:
                    bo = True
                    break
        else:
            bo = True
        if not bo:
            name = [i[0],]
            if len(i) > 5: name += [k for k in i[5:]]
            prolist.append([name,data])
        else:
            raise ValueError('Error: wrong defined entry: {:}'.format(i))
    
    return prolist



def func_calc_M(S):
    """
    Use molecules structure/symbol to calculate molecular weight

    Parameter:
        S  :   structrue in a format: (atomType number) separated by '-' or blank space
               number of '-' and spaces does not matter
               precendent: '-' > blank space 

               Example 1:
                           C2H3O4N5
               Example 2: 
                           C2 - H3 - O4 - N5
               Example 3: 
                           C2   H3   O4   N5
               Example 4: 
                           C2   H3 - O4 - N5
    Return:
        M  : molecular weight (g/mol)
    """
    ##Test list
    ##Slist = [ 123, '  ', '- - ', '---', '1,2,','1 +','4 $',   #bad
    #          'C3H4O5Br1Cl2', 'CHOBrCl','Br Br BrBr',          #good
    #          'C3 - H -2 - 2 - O', 'C3 - H2  2 - O'            #bad]
    log = {'nice':True, }
    # define Periodic Table
    PT = { 'H':1.008, 'B':10.81, 'C':12.01, 'N':14.01, 'O':16.00, 'F':19.00,
           'P':30.91, 'S':32.06, 'Cl':35.45, 'Br':79.90, 'I':126.90 }
    if not isinstance(S,str):
        log['nice'] = False
        log['info'] = 'Error: Molecule structure has to be a string'
        return log, 0.0
    S = S.lower()
    
    proS = []
    # format: split by '-' then split by blank space
    for t in S.split('-'): proS += t.split()
    if len(proS) == 0:
        log['nice'] = False
        log['info'] = 'Error: empty inputs'
        return log, 0.0
    
    proSS = []
    # 1D: split to [  character   number   character   number  ]
    for t in proS:
        if t.isdigit():
            proSS.append(int(t))
        elif t.isalpha():
            proSS.append(t)
        elif t.isalnum():
            stmp = ''
            for c in t:
                if c.isdigit():
                    if stmp.isalpha():
                        proSS.append(stmp)
                        stmp = ''
                else:
                    if stmp.isdigit():
                        proSS.append(int(stmp))
                        stmp = ''
                stmp += c
            if stmp.isdigit():
                proSS.append(int(stmp))
            else:
                proSS.append(stmp)
        else:
            log['nice'] = False
            log['info'] = 'Error: input < {:} > is not correctly defined'.format(t)
            return log, 0.0

    proSSS = []
    # 1D: split to [  atomtype   number   atomtype   number  ]
    for t in proSS:
        if isinstance(t,int):
            proSSS.append(t)
        else:
            # for character, it may have special cases like Br, Cl
            while True:
                if 'br' in t or 'cl' in t:
                    ndx = t.find('br') if 'br' in t else t.find('cl')
                    if ndx > 0: proSSS += [ c for c in t[:ndx] ]
                    proSSS.append(t[ndx:ndx+2])
                    if len(t) >= ndx + 2:
                        t = t[ndx+2:]
                    else:
                        proSSS += [ c for c in t ]
                        break
                else:
                    proSSS += [ c for c in t ]
                    break
    # No adjacent numbers is allowed
    # However the number of each adjacent character is defined at 1
    # Consider cases like:
    #   C 1 2 H <bad> 
    #   C C C 3 <good>
    #   C 1 H 3 <good>
    if not isinstance(proSSS[0],str):
        log['nice'] = False
        log['info'] = 'Error: the atomtype has to be in the first input along with its numbers\n' + \
                      '     : < {:} > is not correctly defined'.format(proSSS[0])
        return log, 0.0
    bo = False
    for t in proSSS:
        if isinstance(t,int):
            if bo:
                log['nice'] = False
                stmp = t
                break
            bo = True
        else:
            bo = False
    if not log['nice']:
        log['info'] = 'Error: no adjacent number inputs is allowd\n' + \
                      '     : < {:} > is not correctly defined'.format(stmp)
        return log, 0.0

    i = 0
    proSSSS = []
    # 2D: [ [atomtype, number],  [atomtype, number], ... ]
    while i < len(proSSS):
        j = i + 1
        if j < len(proSSS) and isinstance(proSSS[j],int):
            proSSSS.append([proSSS[i],proSSS[j]])
            i = j
        else:
            proSSSS.append([proSSS[i],1])
        i += 1
    
    # time to check for Periodic Table
    M = 0.0
    for t in proSSSS:
        tmp = t[0].capitalize()
        if tmp in PT:
            M += PT[tmp] * t[1]
        else:
            log['nice'] = False
            log['info'] = 'Error: atomtype < {:} > is not defined in Periodic Table'.format(tmp)
            break

    return log, M



def func_profile(file,fsize=5):
    """This function is used to process input file
       Return:
            log, fdict"""
    log = {'nice':True,}
    try:
        sizetmp = os.stat(file).st_size
        if sizetmp/1024/1024 > fsize:
            log['nice'] = False
            log['info'] = 'Error: the file size is far larger than {:} MB'.format(fsize)
    except IOError:
        log['nice'] = False
        log['info'] = 'Error : cannot open the file < {:} >!\n'.format(file)
    
    fdict = {}
    if log['nice']:
        profile = []
        # format 2D: [  [ keyword-lower()-strip(), value-strip() ],  ]
        with open(file,'rt') as f:
            while True:
                line = f.readline()
                if len(line) == 0: break
                
                line = line.replace('\n', '')

                proline = line
                if proline.find('#') != -1: proline = proline[:proline.find('#')]
                proline = proline.replace('\t',' ').strip()
                if len(proline) == 0: continue
                
                if proline.find(':') == -1: log['nice'] = False
                if log['nice']:
                    proline = proline.replace('"','').replace("'",'')
                    ltmp = proline.split(':',maxsplit=1)
                    if len(ltmp[0].split()) == 0:
                        log['nice'] = False
                    else:
                        profile.append([ltmp[0].lower().strip(), ltmp[1].strip()])

                if not log['nice']:
                    log['info'] = 'Error in line: < {:} >'.format(line)
                    break

    if log['nice']:
        # get 'mark' indices
        ndxlist = []
        for i,j in enumerate(profile):
            if j[0] == 'mark': ndxlist.append(i)
        if len(ndxlist) > 0 and ndxlist[-1] >= len(profile) - 1:
            log['nice'] = False
            log['info'] = 'Error wrong input file'

        if log['nice']:
            ndxlist.append(len(profile))
            if ndxlist[0] != 0: ndxlist.insert(0,0)

            i = 1
            while i < len(ndxlist):
                beg = ndxlist[i-1]
                end = ndxlist[i]
                i += 1

                key = ''
                tmpdict = {}
                for t in profile[beg:end]:
                    if t[0] in fdict:
                        log['nice'] = False
                        log['info'] = 'Error: double defined < {:} > entry'.format(t[0])
                        break

                    if t[0] == 'mark':
                        if len(t[1]) == 0:
                            log['nice'] = False
                            log['info'] = 'Error: wrong defined < mark > entry < {:} >'.format(t)
                            break
                        key = t[1].lower()
                    elif t[0] == 'type':
                        if len(t[1]) > 0:
                            if t[1].lower() in ['purity', 'pure', 'p']:
                                fdict['type'] = 'purity'
                            elif t[1].lower() in ['mixture', 'mix', 'm']:
                                fdict['type'] = 'mixture'
                            else:
                                log['nice'] = False
                                log['info'] = 'Error: wrong defined < type > entry < {:} >'.format(t)
                                break
                    elif t[0] == 'ratio':
                        if t[1].find(':') == -1:
                            if len(t[1].split()) == 0:
                                y1 = 1.0
                                y2 = 1.0
                            else:
                                log['nice'] = False
                        else:
                            y = t[1].replace(' ','').split(':',maxsplit=1)
                            try:
                                y1 = float(y[0])
                                y2 = float(y[1])
                                if round(y1,5) == 0 or round(y2,5) == 0: raise ValueError
                            except ValueError:
                                log['nice'] = False
                        if not log['nice']:
                            log['info'] = 'Error: wrong defined < ratio > entry < {:} >'.format(t)
                            break
                        y1 = round(y1/(y1+y2),3)
                        y2 = 1.0 - y1
                        fdict['ratio'] = '{:} : {:}'.format(y1,y2)
                        fdict['y1'] = y1
                        fdict['y2'] = y2
                    elif t[0] == 'name':
                        if len(t[1].split()) != 0:
                            name = ''
                            for s in t[1].split(): name += ' ' + s
                            # remember to strip out of blank spaces!!!
                            fdict['name'] = name.strip()
                    elif t[0] in ['m','m1','m2','t-ref','q-ref','t']:
                        if len(t[1]) != 0:
                            try:
                                mw = float(t[1])
                                if mw <= 0: raise ValueError
                            except ValueError:
                                log['nice'] = False
                                log['info'] = 'Error: wrong defined entry < {:} >'.format(t[1])
                                break
                            fdict[t[0]] = mw
                    elif t[0] in ['s','s1','s2']:
                        if len(t[1]) != 0: fdict[t[0]] = t[1]
                    elif t[0] in tmpdict:
                        log['nice'] = False
                        log['info'] = 'Error: double defined group < {:} >'.format(t[0])
                        break
                    else:
                        if len(t[1]) != 0:
                            try:
                                v = int(t[1])
                            except ValueError:
                                log['nice'] = False
                                log['info'] = 'Error: For group < {:} >, < {:} > is not a number'.format(t[0],t[1])
                                break
                        else:
                            v = 1
                        tmpdict[t[0]] = v
                if log['nice']:
                    if len(tmpdict) == 0 and len(key) == 0:
                        pass
                    elif len(tmpdict) != 0 and len(key) != 0:
                        fdict[key] = tmpdict
                    else:
                        log['nice'] = False
                        if len(key) == 0:
                            log['info'] = 'Error: wrong defined entry < {:} >'.format(tmpdict)
                        else:
                            log['info'] = 'Error: no defined groups under < mark > entry'
                        break

    # precendent: type > ring*, no-ring* > m* > s*
    if log['nice']:
        if 'type' not in fdict:
            if 'ring' in fdict or 'no-ring' in fdict:
                fdict['type'] = 'purity'
            elif 'ring-1' in fdict or 'ring-2' in fdict or 'no-ring-1' in fdict or 'no-ring-2' in fdict:
                fdict['type'] = 'mixture'
            else:
                log['nice'] = False
                log['info'] = 'Error: cannot define calculation < type >'

    if log['nice']:
        if fdict['type'] == 'purity':
            if 'm' not in fdict and 's' not in fdict: log['nice'] = False
        if fdict['type'] == 'mixture':
            if 'm1' not in fdict and 's1' not in fdict: log['nice'] = False
            if 'm2' not in fdict and 's2' not in fdict: log['nice'] = False
        if not log['nice']:
            log['info'] = 'Error: the molecular weight is not correctly defined'

    if log['nice']:
        # For future update
        namespace_both = ['name','type','t-ref','q-ref','t']
        namespace_purity = ['ring', 'no-ring', 's', 'm']
        namespace_mixture = ['ring-1','ring-2','no-ring-1','no-ring-2','y1','y2','ratio','s1','s2','m1','m2']
        namespace_purity += namespace_both
        namespace_mixture += namespace_both
        if fdict['type'] == 'purity':
            ctype = 'purity'
            namespace_sel = namespace_purity
        else:
            ctype = 'mixture'
            namespace_sel = namespace_mixture
        for i in fdict:
            if i not in namespace_sel:
                log['nice'] = False
                log['info'] = 'Error: for < {:} >, data type < {:} > is not defined'.format(ctype,i)

    return log, fdict



class MLJR(object):
    """
    Modified Lydersen-Joback-Reid method

    Works only for pure solvent or binary mixture 
    """
    # define constant values
    CONST_AM = 0.5703
    CONST_BM = 1.0121
    CONST_CM = 0.2573
    CONST_EM = 6.75

    def __init__(self, *args, **kwargs):
        self.log = {'nice':True, }
        if 'type' in kwargs:
            self.type = kwargs['type'] 
        else:
            self.log['nice'] = False
            self.log['info'] = 'Error: calculation type is not defined'
            return

        if self.type == 'purity':
            # self.table: 2D: [ [deltaTbM    deltaTM    deltaPM   deltaVM  Number], ... ]
            if 'table' in kwargs:
                self.table = kwargs['table'] 
            else:
                self.log['nice'] = False
                self.log['info'] = 'Error: for < purity >, no input data'
                return
            if 'm' in kwargs:
                self.m = kwargs['m'] 
            else:
                self.log['nice'] = False
                self.log['info'] = 'Error: molecular weight is not defined'
                return
        elif self.type == 'mixture':
            # self.table: 2D: [ [deltaTbM    deltaTM    deltaPM   deltaVM  Number], ... ]
            if 'table1' in kwargs and 'table2' in kwargs:
                self.table1 = kwargs['table1']
                self.table2 = kwargs['table2'] 
            else:
                self.log['nice'] = False
                self.log['info'] = 'Error: for < mixture >, input data is not correctly defined'
                return
            if 'm1' in kwargs and 'm2' in kwargs:
                self.m1 = kwargs['m1'] 
                self.m2 = kwargs['m2'] 
            else:
                self.log['nice'] = False
                self.log['info'] = 'Error: molecular weight is not correctly defined'
                return

            self.y1 = kwargs['y1'] if 'y1' in kwargs else None
            self.y2 = kwargs['y2'] if 'y2' in kwargs else None
            if self.y1 is None and self.y2 is None:
                self.y1 = self.y2 = 0.5
            elif self.y1 is None or self.y2 is None:
                self.log['nice'] = False
                self.log['info'] = 'Error: molecular ratio is not correctly defined'
                return
        else:
            self.log['nice'] = False
            self.log['info'] = 'Error: the input calculation type is wrong'
            return


    def func_calc_critical(self, vTbM, nTbM, vTM, nTM, vPM, nPM, vVM, nVM, M):
        """Return critical properties: Tc, Pc, Vc, Tb"""
        Tb = self.func_calc_Tb(vTbM, nTbM)
        Tc = self.func_calc_Tc(Tb, vTM, nTM)
        Pc = self.func_calc_Pc(M, vPM, nPM)
        Vc = self.func_calc_Vc(vVM, nVM)

        return Tc, Pc, Vc, Tb


    def run(self):
        if self.type == 'purity':
            self.purity()
            self.w = self.func_calc_w(self.Tb, self.Tc, self.Pc)
        else:
            self.mixture()


    def purity(self):
        """
        For pure solvent

        Note: No self-check, if any errors happen, return ValueError Exception

        Return:
            Tc, Pc, Vc, Tb
        """
        # self.table: 2D: [ [deltaTbM    deltaTM    deltaPM   deltaVM   Number], ... ]
        TbM_vlist = [i[0] for i in self.table]
        TbM_nlist = [i[4] for i in self.table]
        TM_vlist = [i[1] for i in self.table]
        TM_nlist = [i[4] for i in self.table]
        PM_vlist = [i[2] for i in self.table]
        PM_nlist = [i[4] for i in self.table]
        VM_vlist = [i[3] for i in self.table]
        VM_nlist = [i[4] for i in self.table]
        self.Tc, self.Pc, self.Vc, self.Tb =  self.func_calc_critical(TbM_vlist, TbM_nlist, TM_vlist, TM_nlist,
                                                                      PM_vlist, PM_nlist, VM_vlist, VM_nlist, self.m)


    def mixture(self):
        """
        For binary mixture solvent

        Note: No self-check, if any errors happen, return ValueError Exception

        Return:
            Tcm, Pcm, Vcm, wm
        """
        # self.table: 2D: [ [deltaTbM    deltaTM    deltaPM   deltaVM   Number], ... ]
        m1_TbM_vlist = [i[0] for i in self.table1]
        m1_TbM_nlist = [i[4] for i in self.table1]
        m1_TM_vlist = [i[1] for i in self.table1]
        m1_TM_nlist = [i[4] for i in self.table1]
        m1_PM_vlist = [i[2] for i in self.table1]
        m1_PM_nlist = [i[4] for i in self.table1]
        m1_VM_vlist = [i[3] for i in self.table1]
        m1_VM_nlist = [i[4] for i in self.table1]

        # get total M
        self.m = self.y1*self.m1 + self.y2*self.m2

        self.Tc1, self.Pc1, self.Vc1, self.Tb1 = self.func_calc_critical(m1_TbM_vlist, m1_TbM_nlist, m1_TM_vlist, m1_TM_nlist,
                                                                         m1_PM_vlist, m1_PM_nlist, m1_VM_vlist, m1_VM_nlist,self.m1)
        self.w1 = self.func_calc_w(self.Tb1, self.Tc1, self.Pc1)

        Pc1_tmp = self.func_calc_Pc(self.m, m1_PM_vlist, m1_PM_nlist)
        w1_tmp = self.func_calc_w(self.Tb1, self.Tc1, Pc1_tmp)

        m2_TbM_vlist = [i[0] for i in self.table2]
        m2_TbM_nlist = [i[4] for i in self.table2]
        m2_TM_vlist = [i[1] for i in self.table2]
        m2_TM_nlist = [i[4] for i in self.table2]
        m2_PM_vlist = [i[2] for i in self.table2]
        m2_PM_nlist = [i[4] for i in self.table2]
        m2_VM_vlist = [i[3] for i in self.table2]
        m2_VM_nlist = [i[4] for i in self.table2]
        self.Tc2, self.Pc2, self.Vc2, self.Tb2 = self.func_calc_critical(m2_TbM_vlist, m2_TbM_nlist, m2_TM_vlist, m2_TM_nlist,
                                                                         m2_PM_vlist, m2_PM_nlist, m2_VM_vlist, m2_VM_nlist,self.m2)
        self.w2 = self.func_calc_w(self.Tb2, self.Tc2, self.Pc2)

        Pc2_tmp = self.func_calc_Pc(self.m, m2_PM_vlist, m2_PM_nlist)
        w2_tmp = self.func_calc_w(self.Tb2, self.Tc2, Pc2_tmp)

        # workaround
        Vcij = self.func_calc_Vcij(self.Vc1,self.Vc2)
        Tcij = self.func_calc_Tcij(self.Tc1,self.Tc2)

        self.Vcm = self.func_calc_Vcm(self.y1,self.y2,Vcij)
        self.Tcm = self.func_calc_Tcm(self.Vcm,self.y1,self.y2,Vcij,Tcij)

        self.wm = self.func_calc_wm(w1_tmp, self.y1, w2_tmp, self.y2)
        self.Pcm = self.func_calc_Pcm(self.Tcm, self.Vcm, self.wm)
        self.Tbm = self.Tb1 * self.y1 + self.Tb2 * self.y2


    def func_calc_sum(self, vlist, nlist=None):
        """
        Parameters:
            vlist: 1D list, float, values to be calculated
            nlist: 1D list, integer, number of each indices

            Note: they have to correspond to each other
        Return:
            SUM{ ni * vi }
        """
        if nlist is None: nlist = [1 for i in range(len(vlist))]

        return sum([v*nlist[i] for i,v in enumerate(vlist)])


    def func_calc_Tb(self, vlist, nlist=None):
        """
        Parameters:
            vlist: 1D list, float, values to be calculated
            nlist: 1D list, float, number of each indices

            Note: they have to correspond to each other
        Return:
            Tb = 198.2 + SUM{ n*deltaTbM }
        """
        return 198.2 + self.func_calc_sum(vlist,nlist)


    def func_calc_Tc(self, Tb, vlist, nlist=None):
        """
        Parameters:
            Tb   : float

            vlist: For TM, 1D list, float, values to be calculated
            nlist: For TM, 1D list, float, number of each indices

            Note: they have to correspond to each other
        Return:
            Tc = Tb / [ AM + BM * SUM(n*deltaTM) - (SUM(N*deltaTM))^2 ]
        """
        #calc SUM(n*deltaTM)
        STM = self.func_calc_sum(vlist,nlist)

        return Tb / (self.CONST_AM + self.CONST_BM * STM - STM * STM)


    def func_calc_Pc(self, M, vlist, nlist=None):
        """
        Parameters:
            vlist: 1D list, float, values to be calculated
            nlist: 1D list, float, number of each indices

            M    : molecular weight (unit in g/mol)

            Note: they have to correspond to each other
        Return:
            Tb = M / ( CM + SUM{ ni * vi } )^2
        """
        SP = self.func_calc_sum(vlist,nlist)
        return M / (self.CONST_CM + SP) / (self.CONST_CM + SP)


    def func_calc_Vc(self, vlist, nlist=None):
        """
        Parameters:
            vlist: 1D list, float, values to be calculated
            nlist: 1D list, float, number of each indices

            Note: they have to correspond to each other
        Return:
            Tb = EM + SUM{ ni * vi }
        """
        return self.CONST_EM + self.func_calc_sum(vlist,nlist)


    def func_calc_Vcm(self, yi, yj, Vcij):
        """
        Parameters:
            yi   : float, ratio for component i  
            yj   : float, ratio for component j
            Vcij : 2D matrix, in C-type, [ [Vcii, Vcij], [ Vcji, Vcjj] ]

            Note: they have to correspond to each other
        Return:
            Vcm = SUM{SUM{yi*yj*Vcij}}
        """
        yii = yi*yi*Vcij[0][0]
        yij = yi*yj*Vcij[0][1]
        yji = yj*yi*Vcij[1][0]
        yjj = yj*yj*Vcij[1][1]

        return yii + yij + yji + yjj


    def func_calc_Vcij(self, Vci, Vcj):
        """
        Parameters:
            Vci : float, for component i
            Vcj : float, for component j

            Note: they have to correspond to each other
        Return:
            Vcij : 2D matrix, in C-type, [ [Vcii, Vcij], [ Vcji, Vcjj] ]
        """
        def func_calc_cube(i,j):
            return pow(pow(i,1/3)+pow(j,1/3),3)
        
        return [[1/8*func_calc_cube(Vci,Vci), 1/8*func_calc_cube(Vci,Vcj)],
                [1/8*func_calc_cube(Vcj,Vci), 1/8*func_calc_cube(Vcj,Vcj)]]


    def func_calc_Tcij(self, Tci, Tcj):
        """
        Parameters:
            Tci : float, for component i
            Tcj : float, for component j

            Note: kij is assumed to be 1 due to lack of literature values
            Note: they have to correspond to each other
        Return:
            Tcij : 2D matrix, in C-type, [ [Tcii, Tcij], [ Tcji, Tcjj] ]
        """
        def func_calc_power(i,j):
            return pow(i*j,1/2)
        
        return [[func_calc_power(Tci,Tci), func_calc_power(Tci,Tcj)],
                [func_calc_power(Tcj,Tci), func_calc_power(Tcj,Tcj)]]


    def func_calc_Tcm(self, Vcm, yi, yj, Vcij, Tcij):
        """
        Parameters:
            Vcm  : float, for binary mixture, critical Volume
            yi   : float, ratio for component i  
            yj   : float, ratio for component j
            Vcij : 2D matrix, in C-type, [ [Vcii, Vcij], [ Vcji, Vcjj] ]
            Tcij : 2D matrix, in C-type, [ [Tcii, Tcij], [ Tcji, Tcjj] ]

            Note: they have to correspond to each other
        Return:
            Tcm = 1 / pow(Vcm,1/4) * SUM{SUM{yi*yj*pow(Vij,1/4)*Tcij}}
        """
        tii = yi * yi * pow(Vcij[0][0],1/4) * Tcij[0][0]
        tij = yi * yj * pow(Vcij[0][1],1/4) * Tcij[0][1]
        tji = yj * yi * pow(Vcij[1][0],1/4) * Tcij[1][0]
        tjj = yj * yj * pow(Vcij[1][1],1/4) * Tcij[1][1]

        return (tii + tij + tji + tjj) / pow(Vcm,1/4)


    def func_calc_w(self, Tb, Tc, Pc, Pb=1.0):
        """
        Parameters:
            Tb  : float, for pure solvent, normal boiling temperature
            Tc  : float, for pure solvent, critical temperature
            Pc  : float, for pure solvent, critical pressure

            Pb  : 1 bar (from literature testing, assume 1 atm = 1 bar = 100000 Pa)
        Return:
            w: arentric factor for pure solvent
        """
        t = (Tb-43)*(Tc-43) / ((Tc-Tb)*(0.7*Tc-43))
        l = math.log10(Pc/Pb)
        s = (Tc-43) / (Tc-Tb)

        return t*l - s*l + l - 1


    def func_calc_wm(self, wi, yi, wj, yj):
        """
        Parameters:
            wi  : float, for binary mixture, arentric factor for i
            yi  : float, ratio for component i

            wj  : float, for binary mixture, arentric factor for j
            yj  : float, ratio for component j
        Return:
            wm: arentric factor for binary mixture
        """
        return wi*yi + wj*yj


    def func_calc_Pcm(self, Tcm, Vcm, wm, R=8.314):
        """
        Parameters:
            Tcm  : float, for binary mixture, critical Temperature
            Vcm  : float, for binary mixture, critical Volume
            wm   : float, for binary mixture, acentric factor

            R    : Gas constant, 8.314 J/mol K
        Return:
            Paper1 has the error, correct is:
            Pcm = (0.2905 - 0.085*wm) * R * Tcm / Vcm
        """
        return ( 0.2905 - 0.085*wm ) * 8.314 * Tcm / Vcm * 10



def func_calc_st(T,Qref,Tref,Tc):
    """
    Calculate Elevated Temperature Surface Tension
    Parameters:
        Qref  :  Reference Surface Tension at temperature T, (mN/m)
        Tref  :  Reference temperature (K)
        Tc    :  Critical temperature (K)
    Return:
        Surface Tension at temperature T
    """
    tmp = (Tc-T)/(Tc-Tref)
    return pow(tmp,11/9)*Qref


def func_calc_density(M,T,Tc,Pc,Vc,Tb):
    """
    Function to calculate density

    Parameters:
        Tc  :  Critical temperature (K)
        Pc  :  Critical pressure (bar)
        Vc  :  Critical molar volume (cm^3/mol)
        Tb  :  Normal boiling temperature (K)
        M   :  Molecular weight (g/mol)
        
    Return:
        Be extremely careful about the unit: 1 bar = 100000 Pa
        D   :  density (g/mL)
    """
    Tr = T / Tc
    Tbr = Tb / Tc
    beta = (-1.0 - pow(1-Tr,2/7)) / (1 + pow(1-Tbr, 2/7))
    tmp = pow(0.3445*Pc*pow(Vc,1.0135)/Tc/8.314/10, beta)
    D = M * Pc / 8.314 / Tc * tmp / 10
    return D



def func_pro_argparse():
    """
    Process input arguments
    Return:
       Namespace(file=None, sample=False, template=False, CCG=False, CCl=False, avail_group_name=False, output=False)
    """
    parser = argparse.ArgumentParser(description='Critical Properties Calculation using Modified Lydersen-Joback-Reid method',allow_abbrev=False)
    parser.add_argument('-v','--version',action='version',version='mLJR {:}'.format(__version__))
    parser.add_argument('-f','--file',help='Input file path')
    parser.add_argument('-s','--sample',help='Template shows on screen (default)',action='store_true')
    parser.add_argument('-t','--template',help='Template writes into a file',action='store_true')
    parser.add_argument('--CCl',help='Template molecule, Choline Chloride <purity> (default)',action='store_true')
    parser.add_argument('--CCG',help='Template molecule, Choline Chloride Glycerol <mixture>',action='store_true')
    parser.add_argument('-g','--avail-group-name',help='Available Group Names <include nick names>',action='store_true')
    parser.add_argument('-o','--output',help='Append results to input file (default False)',action='store_true')
    parser.add_argument('-e','--examples',help='Show command line examples',action='store_true')
    parser.add_argument('-x','--calc',help='If only one input, it will be thought as T, the density value will be \
                            calculated, if it has three inputs, which has to be input in a sequence (T-ref, Q-ref, T), \
                            then both density and surface tension will be calculated, if else, wrong',
                            nargs='+',type=float)

    return parser



def main():
    """main!"""

    # process arguments
    parser = func_pro_argparse()
    args = parser.parse_args()
    #Namespace(file=None, sample=False, template=False, CCG=False, CCl=False, avail_group_name=False, output=False)
    # avoid no inputs, print help message instead
    bo = True
    for i in args.__dict__:
        if not ( args.__dict__[i] is False or args.__dict__[i] is None ):
            bo = False
            break
    if bo:
        parser.print_help()
        exit()
    
    # particular check density and surface tension inputs
    if args.calc is not None and len(args.calc) not in [1,3]:
        print('Error: the input parameter(s) in -x/--calc is not correctly defined')
        exit()


    if args.examples:
        txt = '# Python version 3\n'
        txt += '# If you have `python3` in your command, directly use:\n'
        txt += '#    mljr  [options]\n'
        txt += '# If not, you have to run this script in:\n'
        txt += '#    [python3] mljr  [optioins]\n\n'
        txt += '# To check all avail_group_name:\n'
        txt += '#    [python3] mljr  -g\n\n'
        txt += '# Show default template to screen:\n'
        txt += '#    [python3] mljr  -s\n\n'
        txt += '# Generate default template (CCl) to a file:\n'
        txt += '#    [python3] mljr  -t\n\n'
        txt += '# Or do Both\n'
        txt += '#    [python3] mljr  -s  -t\n\n'
        txt += '# Similarly, for another template, show it on screen:\n'
        txt += '# For -s -t --CCG --CCl, you can combine them as many as you want\n'
        txt += '#    [python3] mljr  -s  --CCG\n\n'
        txt += '# For calculation, show results to screen (default)\n'
        txt += '#    [python3] mljr  -f [file]\n\n'
        txt += '# For calculation, append results to the input file\n'
        txt += '#    [python3] mljr  -f [file] -o\n\n'
        print(txt)
        exit()


    pro_data_ljr_ring = func_proDATA(DATA_LJR_RING)
    pro_data_ljr_noring = func_proDATA(DATA_LJR_NO_RING)
    # if chosen, print available group names
    if args.file is None:
        if args.avail_group_name:
            print('Quoted in brace, may have many names, case-insensitive')
            print('Symbol @ means three bonds\n')
            print('For Ring, Available Group Names are:')
            tot = '--> '
            for i in pro_data_ljr_ring:
                line = ''
                for j in i[0]: line += j + '   '
                line = '{ ' + line.strip() + ' }'

                tot += line + '    '
                if len(tot) > 54:
                    print(tot.strip())
                    tot = '--> '
            print(tot)

            print('\nFor No-Ring, Available Group Names are:')
            tot = '--> '
            for i in pro_data_ljr_noring:
                line = ''
                for j in i[0]: line += j + '   '
                line = '{ ' + line.strip() + ' }'

                tot += line + '    '
                if len(tot) > 54:
                    print(tot.strip())
                    tot = '--> '
            print(tot)

        if args.sample or args.template:
            if not args.CCl and not args.CCG: args.CCl = True
        
        if args.CCl or args.CCG:
            if not args.sample and not args.template: args.sample = True

        if args.sample:
            if args.CCl: func_show_template(TEMPLATE_CCL)
            if args.CCG: func_show_template(TEMPLATE_CCG)
        
        if args.template:
            if args.CCl: func_show_template(TEMPLATE_CCL,file=True)
            if args.CCG: func_show_template(TEMPLATE_CCG,file=True)
        
        exit()
    elif not os.path.isfile(args.file):
        print('Error: wrong input file < {:} >'.format(args.file))
        exit()
    
    # process input file, at the same time get their molecular weight
    log, fdict = func_profile(args.file)
    if not log['nice']:
        print(log['info'])
        exit()

    if fdict['type'] == 'purity':
        if 'm' not in fdict:
            log, m = func_calc_M(fdict['s'])
            if log['nice']:
                fdict['m'] = m
            else:
                print(log['info'])
                exit()
        table = []
        if 'ring' in fdict:
            for sym, nm in fdict['ring'].items():
                bo = True
                for ref in pro_data_ljr_ring:
                    if sym in [i.lower() for i in ref[0]]:
                        # Note: make a copy!!
                        table.append( ref[1][:] + [nm,] )
                        bo = False
                        break
                if bo:
                    raise ValueError('Error: for ring, symbol < {:} > is not defined'.format(sym))
        if 'no-ring' in fdict:
            for sym, nm in fdict['no-ring'].items():
                bo = True
                for ref in pro_data_ljr_noring:
                    if sym in [i.lower() for i in ref[0]]:
                        # Note: make a copy!!
                        table.append( ref[1][:] + [nm,] )
                        bo = False
                        break
                if bo:
                    raise ValueError('Error: for no-ring, symbol < {:} > is not defined'.format(sym))
        fdict['table'] = table

    else:
        if 'm1' not in fdict:
            log, m1 = func_calc_M(fdict['s1'])
            if log['nice']:
                fdict['m1'] = m1
            else:
                print(log['info'])
                exit()
        if 'm2' not in fdict:
            log, m2 = func_calc_M(fdict['s2'])
            if log['nice']:
                fdict['m2'] = m2
            else:
                print(log['info'])
                exit()

        #[    [ [ allNames ], [ deltaTbM    deltaTM    deltaPM   deltaVM ] ],    ]
        table1 = []
        if 'ring-1' in fdict:
            for sym, nm in fdict['ring-1'].items():
                bo = True
                for ref in pro_data_ljr_ring:
                    if sym in [i.lower() for i in ref[0]]:
                        # Note: make a copy!!
                        table1.append( ref[1][:] + [nm,] )
                        bo = False
                        break
                if bo:
                    raise ValueError('Error: for ring, symbol < {:} > is not defined'.format(sym))
        if 'no-ring-1' in fdict:
            for sym, nm in fdict['no-ring-1'].items():
                bo = True
                for ref in pro_data_ljr_noring:
                    if sym in [i.lower() for i in ref[0]]:
                        # Note: make a copy!!
                        table1.append( ref[1][:] + [nm,] )
                        bo = False
                        break
                if bo:
                    raise ValueError('Error: for no-ring, symbol < {:} > is not defined'.format(sym))
        table2 = []
        if 'ring-2' in fdict:
            for sym, nm in fdict['ring-2'].items():
                bo = True
                for ref in pro_data_ljr_ring:
                    if sym in [i.lower() for i in ref[0]]:
                        # Note: make a copy!!
                        table2.append( ref[1][:] + [nm,] )
                        bo = False
                        break
                if bo:
                    raise ValueError('Error: for ring, symbol < {:} > is not defined'.format(sym))
        if 'no-ring-2' in fdict:
            for sym, nm in fdict['no-ring-2'].items():
                bo = True
                for ref in pro_data_ljr_noring:
                    if sym in [i.lower() for i in ref[0]]:
                        # Note: make a copy!!
                        table2.append( ref[1][:] + [nm,] )
                        bo = False
                        break
                if bo:
                    raise ValueError('Error: for no-ring, symbol < {:} > is not defined'.format(sym))
        fdict['table1'] = table1
        fdict['table2'] = table2


    rst = MLJR(**fdict)

    if rst.log['nice']:
        rst.run()
    else:
        print(rst.log['info'])
        exit()

    bo_calc_st = True
    bo_calc_dy = True
    if args.calc is not None:
        if len(args.calc) == 1:
            t = args.calc[0]
            bo_calc_st = False
        else:
            t_ref = args.calc[0]
            q_ref = args.calc[1]
            t = args.calc[2]
    elif 't-ref' in fdict:
        t_ref = fdict['t-ref']
        if 'q-ref' in fdict and 't' in fdict:
            q_ref = fdict['q-ref']
            t = fdict['t']
        else:
            bo_calc_st = False
    else:
        bo_calc_st = False
        bo_calc_dy = False

    info = '\n\n'
    if rst.type == 'purity':
        info += '# For calculation type < purity >\n'
        info += '# Molecular weight                   m = < {:} >\n\n'.format(round(rst.m,4))
        info += '# Critical temperature (K):          Tc = < {:} >\n'.format(round(rst.Tc,4))
        info += '# Critical pressure (bar):           Pc = < {:} >\n'.format(round(rst.Pc,4))
        info += '# Critical molar volume (cm^3/mol):  Vc = < {:} >\n'.format(round(rst.Vc,4))
        info += '# Boiling temperature (K):           Tb = < {:} >\n'.format(round(rst.Tb,4))
        info += '# Acentric factor:                   w  = < {:} >\n\n'.format(round(rst.w,4))
        if bo_calc_dy:
            den = func_calc_density(rst.m,t,rst.Tc,rst.Pc,rst.Vc,rst.Tb)
            info += '# Density at ({:}):           d  = < {:} >\n\n'.format(round(t,2),round(den,4))
        if bo_calc_st:
            st = func_calc_st(t,q_ref,t_ref,rst.Tc)
            info += '# Surf. Ten. at ({:}):       st  = < {:} >\n\n'.format(round(t,2),round(st,4))
    else:
        info += '# For calculation type < mixture >\n'
        info += '# Molar ratio: ( m1 : m2 ) = < {:} : {:} >\n'.format(round(rst.y1,4),round(rst.y2,4))
        info += '# Total Molecular weight:   m = < {:} >\n\n\n'.format(round(rst.m,4))

        info += '# Molecular weight                 m1  = < {:} >\n'.format(round(rst.m1,4))
        info += '# Molar ratio:                     y1  = < {:} >\n'.format(round(rst.y1,4))
        info += '# Critical temperature (K):        Tc1 = < {:} >\n'.format(round(rst.Tc1,4))
        info += '# Critical pressure (bar):         Pc1 = < {:} >\n'.format(round(rst.Pc1,4))
        info += '# Critical molar volume (mL/mol):  Vc1 = < {:} >\n'.format(round(rst.Vc1,4))
        info += '# Boiling temperature (K):         Tb1 = < {:} >\n'.format(round(rst.Tb1,4))
        info += '# Acentric factor:                 w1  = < {:} >\n\n'.format(round(rst.w1,4))

        info += '# Molecular weight                 m2  = < {:} >\n'.format(round(rst.m2,4))
        info += '# Molar ratio:                     y2  = < {:} >\n'.format(round(rst.y2,4))
        info += '# Critical temperature (K):        Tc2 = < {:} >\n'.format(round(rst.Tc2,4))
        info += '# Critical pressure (bar):         Pc2 = < {:} >\n'.format(round(rst.Pc2,4))
        info += '# Critical molar volume (mL/mol):  Vc2 = < {:} >\n'.format(round(rst.Vc2,4))
        info += '# Boiling temperature (K):         Tb2 = < {:} >\n'.format(round(rst.Tb2,4))
        info += '# Acentric factor:                 w2  = < {:} >\n\n\n'.format(round(rst.w2,4))

        info += '# Mixing Critical temperature (K):        Tcm = < {:} >\n'.format(round(rst.Tcm,4))
        info += '# Mixing Critical pressure (bar):         Pcm = < {:} >\n'.format(round(rst.Pcm,4))
        info += '# Mixing Critical molar volume (mL/mol):  Vcm = < {:} >\n'.format(round(rst.Vcm,4))
        info += '# Mixing Boiling temperature (K):         Tbm = < {:} >\n'.format(round(rst.Tbm,4))
        info += '# Mixing Acentric factor:                 wm  = < {:} >\n\n'.format(round(rst.wm,4))

        if bo_calc_dy:
            den = func_calc_density(rst.m,t,rst.Tcm,rst.Pcm,rst.Vcm,rst.Tbm)
            info += '# Density   at  {:} K (g/mL):    d = < {:} >\n'.format(round(t,2),round(den,4))
        if bo_calc_st:
            st = func_calc_st(t,q_ref,t_ref,rst.Tcm)
            info += '# Surf. Ten. at {:} K (mN/m):   st = < {:} >\n'.format(round(t,2),round(st,4))
            info += '# Reference: Temp. < {:} >, Surface Tension: < {:} >\n\n'.format(round(t_ref,4),round(q_ref,4))
        

    if args.output:
        with open(args.file, 'a+') as f: f.write(info)
    else:
        print(info,end='')

if __name__ == '__main__':
    main()


