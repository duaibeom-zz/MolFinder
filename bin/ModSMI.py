import os
import re
import time
import random
import itertools

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

limit_ = 200


def set_avoid_ring(_smiles):
    avoid_ring = []
    ring_tmp = set(re.findall(r"\d", _smiles))
    for j in ring_tmp:
        tmp = [i for i, val in enumerate(_smiles) if val == j]
        while tmp:
            avoid_ring += [j for j in range(tmp.pop(0), tmp.pop(0) + 1)]
    return set(avoid_ring)


def prepare_rigid_crossover(_smiles, side, consider_ring=True, _minimum_len=4):
    """ 1 point crossover
    :param _smiles: SMILES (str)
    :param side: Left SMILES or Right SMILES ['L'|'R'] (str)
    :param consider_ring: consider about avoiding ring (bool)
    :param _minimum_len: minimum cut size (int)
    :return:
    """

    if side not in ["L", "R"]:
        raise Exception("You must choice in L(Left) or R(Right)")

    _smiles_len = len(_smiles)
    _smi = None

    if consider_ring:  # ring 부분을 피해서 자르도록 고려.
        avoid_ring_list = set_avoid_ring(_smiles)

    p = 0
    _start = None
    _end = None
    _gate = False
    while not _gate:  # 원하는 형태의 smiles piece가 나올 때 까지 반복
        if p == limit_:
            raise ValueError(f"main_gate fail ({side}): {_smiles}")

        if consider_ring:
            j = 0
            ring_gate = False
            if side == "L":
                while not ring_gate:
                    if j == 30:
                        raise ValueError(f"ring_gate fail (L): {_smiles}")
                    _end = np.random.randint(_minimum_len, _smiles_len + 1)
                    if _end not in avoid_ring_list:
                        ring_gate = True
                    j += 1
            elif side == "R":
                while not ring_gate:
                    if j == 30:
                        raise ValueError(f"ring_gate fail (R): {_smiles}")
                    _start = np.random.randint(0, _smiles_len - _minimum_len)
                    if _start not in avoid_ring_list:
                        ring_gate = True
                    j += 1
            _smi = _smiles[_start:_end]
        else:
            if side == "L":
                _end = np.random.randint(_minimum_len, _smiles_len)
            elif side == "R":
                _start = np.random.randint(0, _smiles_len - _minimum_len)

            _smi = _smiles[_start:_end]
            chk_ring = re.findall(r"\d", _smi)
            i = 0
            for i in set(chk_ring):
                list_ring = [_ for _, val in enumerate(_smi) if val == i]
                if (len(list_ring) % 2) == 1:
                    b = random.sample(list_ring, 1)
                    _smi = _smi[:b[0]] + _smi[b[0] + 1:]
                    # print(f'@ {_smi} // {_smiles}')

        p += 1

        if "." in _smi:  # 이온은 패스한다.
            continue

        n_chk = 0
        for j in _smi:  # [] 닫혀 있는지 확인.
            if j == "[":
                n_chk += 1
            if j == "]":
                n_chk -= 1
        if n_chk == 0:
            _gate = True

    return _smi


def chk_branch(_smi, side=None):

    if side not in ["L", "R", None]:
        raise Exception("You must choice in L(Left) or R(Right)")

    branch_list = []
    n_branch = 0
    min_branch = 0
    for i, b in enumerate(_smi):  # () 닫혀 있는지 확인.
        if b == "(":
            n_branch += 1  # 0 == (
        if b in ")":
            n_branch -= 1  # 1 == )
        if side == "L":
            if n_branch < min_branch:
                min_branch = n_branch
                branch_list.append(i)
        elif side == "R":  # track max_value
            if n_branch > min_branch:
                min_branch = n_branch
                branch_list.append(i)
    if side == None:
        return n_branch
    return np.asarray(branch_list), min_branch


def get_open_branch(_smi):
    return [i for i, e in enumerate(_smi) if e == "("]


def get_close_branch(_smi):
    return [i for i, e in enumerate(_smi) if e == ")"]


def tight_rm_branch(_smi_l, _smi_r):
    # tmp = time.time()

    _new_smi = _smi_l + _smi_r

    open_branch = get_open_branch(_new_smi)
    close_branch = get_close_branch(_new_smi)

    b = None
    n_branch = chk_branch(_new_smi)

    q = len(_smi_l)
    while n_branch > 0:  # over opened-branch
        _smi_l_open_branch = get_open_branch(_smi_l)
        _smi_r_open_branch = get_open_branch(_smi_r)
        open_branch = get_open_branch(_smi_l + _smi_r)
        avoid_tokens = [
            i for i, e in enumerate(_smi_l + _smi_r)
            if e in ["=", "#", "@", "1", "2", "3", "4", "5", "6", "7", "8"]
        ]

        if len(_smi_r_open_branch) == 0:  # open branch 가 없을 경우
            _smi_r_open_branch.append(len(_smi_r))
        if len(_smi_l_open_branch) == 0:
            _smi_l_open_branch.append(0)

        n = np.random.rand()  # 임의적으로 close branch 를 추가하거나 제거한다.
        if n > 0.5:  # 추가
            branch_gate = False
            j = 0
            while not branch_gate:  # Ring 부분을 피해서 자름
                if j == limit_:
                    raise ValueError
                b = np.random.randint(_smi_l_open_branch[-1] + 1,
                                      _smi_r_open_branch[-1] + q)
                j += 1
                if b not in avoid_tokens:
                    branch_gate = True
            n_branch -= 1
            if b <= len(_smi_l
                        ):  # SMILES 길이를 고려하여 자른다. 좌측 SMILES의 open branch를 cut!
                _smi_l = _smi_l[:b] + ")" + _smi_l[b:]
                q += 1
            else:  # 좌측 SMILES 길이를 제외한 수가 우측 SMILES 문자의 위치를 의미한다.
                b -= len(_smi_l)
                _smi_r = _smi_r[:b] + ")" + _smi_r[b:]
        else:  # 제거
            b = _smi_l_open_branch[-1]  # (Random으로도 가능함. 과한 부분만 Cut!)
            n_branch -= 1
            q -= 1
            _smi_l = _smi_l[:b] + _smi_l[b + 1:]

    while n_branch < 0:  # over closed-branch
        _smi_l_close_branch = get_close_branch(_smi_l)
        _smi_r_close_branch = get_close_branch(_smi_r)
        close_branch = get_close_branch(_smi_l + _smi_r)
        avoid_tokens = [
            i for i, e in enumerate(_smi_l + _smi_r)
            if e in ["=", "#", "@", "1", "2", "3", "4", "5", "6", "7", "8"]
        ]

        if len(_smi_r_close_branch) == 0:
            _smi_r_close_branch.append(len(_smi_r))
        if len(_smi_l_close_branch) == 0:
            _smi_l_close_branch.append(0)

        n = np.random.rand()
        if n > 0.5:
            branch_gate = False
            j = 0
            while not branch_gate:  # Ring 부분을 피해서 자름
                b = np.random.randint(_smi_l_close_branch[-1] + 1,
                                      _smi_r_close_branch[0] + q + 1)
                j += 1
                if b not in (close_branch + avoid_tokens):
                    branch_gate = True
                if j == limit_:
                    raise ValueError
            n_branch += 1
            if b < len(_smi_l):
                _smi_l = _smi_l[:b] + "(" + _smi_l[b:]
                q += 1
            else:
                b -= len(_smi_l)
                _smi_r = _smi_r[:b] + "(" + _smi_r[b:]
        else:
            b = _smi_r_close_branch[0]
            n_branch += 1
            # print(f'{_smi_r[b]}')
            _smi_r = _smi_r[:b] + _smi_r[b + 1:]

    # time_.append(time.time() - tmp)

    return _smi_l + _smi_r


def replace_atom(_smi):

    #                    C /B  N  P / O  S / F  Cl  Br  I
    replace_atom_list = [6, 5, 7, 15, 8, 16, 9, 17, 35, 53]
    #                         C  N  P / O  S
    replace_arom_atom_list = [6, 7, 15, 8, 16]

    # print(f"before: {_smi}")

    mol_ = Chem.MolFromSmiles(_smi)
    max_len = mol_.GetNumAtoms()

    mw = Chem.RWMol(mol_)
    # Chem.SanitizeMol(mw)

    p = 0
    gate_ = False
    while not gate_:
        if p == 30:
            # raise Exception
            raise PermissionError

        rnd_atom = np.random.randint(0, max_len)

        valence = mw.GetAtomWithIdx(rnd_atom).GetExplicitValence()
        if mw.GetAtomWithIdx(rnd_atom).GetIsAromatic():
            if valence == 3:
                mw.ReplaceAtom(
                    rnd_atom,
                    Chem.Atom(replace_arom_atom_list[np.random.randint(0, 3)]))
            elif valence == 2:
                mw.ReplaceAtom(
                    rnd_atom,
                    Chem.Atom(replace_arom_atom_list[np.random.randint(1, 5)]))
            else:
                continue
            mw.GetAtomWithIdx(rnd_atom).SetIsAromatic(True)
        else:
            if valence == 4:
                mw.ReplaceAtom(
                    rnd_atom,
                    Chem.Atom(replace_atom_list[np.random.randint(0, 1)]))
            elif valence == 3:
                mw.ReplaceAtom(
                    rnd_atom,
                    Chem.Atom(replace_atom_list[np.random.randint(0, 4)]))
            elif valence == 2:
                mw.ReplaceAtom(
                    rnd_atom,
                    Chem.Atom(replace_atom_list[np.random.randint(0, 6)]))
            elif valence == 1:
                mw.ReplaceAtom(
                    rnd_atom,
                    Chem.Atom(replace_atom_list[np.random.randint(0, 10)]))

        p += 1
        # print(f"after: {Chem.MolToSmiles(mw)}")
        try:
            Chem.SanitizeMol(mw)
            gate_ = True
        except Chem.rdchem.KekulizeException:
            # print(f"{_smi} {Chem.MolToSmiles(mw, kekuleSmiles=False)} {Chem.MolToSmiles(mol_, kekuleSmiles=False)}")
            # raise Exception
            pass

    Chem.Kekulize(mw)
    # print(f"after_San: {Chem.MolToSmiles(mw)}")

    return Chem.MolToSmiles(mw, kekuleSmiles=True), mw


def delete_atom(_smi):
    """
    Aromatic ring 을 제외하고 삭제함.
    :param _smi:
    :return:
    """
    max_len = len(_smi)
    mol_ = False
    _new_smi = None

    p = 0
    while not mol_:
        p += 1
        if p == 30:
            # raise Exception
            raise PermissionError

        rnd_insert = np.random.randint(max_len)
        _new_smi = _smi[:rnd_insert] + _smi[rnd_insert + 1:]
        mol_ = Chem.MolFromSmiles(_new_smi)

    return _new_smi, mol_


def add_atom(_smi):
    list_atom = ["C", "B", "N", "P", "O", "S", "Cl", "Br"]

    max_len = len(_smi)
    mol_ = False
    _new_smi = None

    p = 0
    while not mol_:
        p += 1
        if p == 30:
            # raise Exception
            raise PermissionError

        rnd_insert = np.random.randint(max_len)
        _new_smi = _smi[:rnd_insert] + random.sample(list_atom,
                                                     1)[0] + _smi[rnd_insert:]
        mol_ = Chem.MolFromSmiles(_new_smi)

    return _new_smi, mol_


def cut_smi(smi1, smi2, func, ring_bool):

    l_smi = None
    r_smi = None

    try:
        l_smi = func(smi1, "L", ring_bool, 4)
        r_smi = func(smi2, "R", ring_bool, 4)
    except (IndexError, ValueError):
        # fail_f.write(f"{l_smi},{r_smi},piece\n")
        raise PermissionError

    return l_smi, r_smi


def crossover_smiles(smi1, smi2, func, ring_bool):
    new_smi = None
    mol = None
    l_smi = None
    r_smi = None

    gate = 0
    while not (l_smi and r_smi):
        gate += 1
        if gate > 10:
            # fail_f.write(f"{l_smi},{r_smi},np\n")
            raise PermissionError
        try:
            l_smi, r_smi = cut_smi(smi1, smi2, func, ring_bool)
        except:
            pass

    gate = 0
    while not mol:
        gate += 1
        if gate > 5:
            break
        try:
            new_smi = tight_rm_branch(l_smi, r_smi)
        except ValueError:
            continue
        mol = Chem.MolFromSmiles(new_smi)

    return new_smi, mol


if __name__ == "__main__":

    from multiprocessing import Pool
    import datetime
    from rkdit import rdBase

    rdBase.DisableLog('rdApp.*')

    _minimum_len = 4

    def get_smi(_smi):
        mol = Chem.MolFromSmiles(_smi[:-1])
        Chem.Kekulize(mol)
        return Chem.MolToSmiles(mol, kekuleSmiles=True, isomericSmiles=False)

    with open(
            '/home/yongbeom/Git/duaibeom/gen_smiles/data/guacamol_v1_test.smiles'
    ) as f:
        with Pool(8) as pool:
            smiles = pool.map(get_smi, f)

    # smiles = smiles[:1000]
    with open('init.smi', 'w') as f:
        for smi in smiles:
            f.write(f'{smi}\n')

    def just_cut_left_smi(_smi):
        len_smi = len(_smi)
        _end = np.random.randint(_minimum_len, len_smi + 1)
        return _smi[:_end]

    def just_cut_right_smi(_smi):
        len_smi = len(_smi)
        _start = np.random.randint(0, len_smi - _minimum_len)
        return _smi[_start:]

    def simple_cut_smi(smi1, smi2):

        l_smi = None
        r_smi = None

        l_smi = just_cut_left_smi(smi1)
        r_smi = just_cut_right_smi(smi2)

        return l_smi + r_smi

    start_time = time.time()
    print("-- start -- simple cut with swap")
    n_valid = 0
    n_fail = 0

    smi_list = []
    for _ in range(100000):
        smi1, smi2 = random.sample(smiles, 2)
        gate = 0
        mol = None
        _smi = None

        while not mol:
            gate += 1
            if gate > 5:
                break
            try:
                _smi = simple_cut_smi(smi1, smi2)
                mol = Chem.MolFromSmiles(_smi)
            except:
                mol = None

        if not mol:
            gate = 0
            while not mol:
                gate += 1
                if gate > 5:
                    break
                try:
                    _smi = simple_cut_smi(smi2, smi1)
                    mol = Chem.MolFromSmiles(_smi)
                except:
                    mol = None

        if mol:
            smi_list.append(_smi)
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")
    # with open('crossover_nothing_swap.smi', 'w') as f:
    #     for smi in smi_list:
    #         f.write(f'{smi}\n')

    start_time = time.time()
    print("-- start -- simple cut")
    n_valid = 0
    n_fail = 0

    smi_list = []
    for _ in range(100000):
        smi1, smi2 = random.sample(smiles, 2)
        gate = 0
        mol = None
        _smi = None

        while not mol:
            gate += 1
            if gate > 5:
                break
            try:
                _smi = simple_cut_smi(smi1, smi2)
                mol = Chem.MolFromSmiles(_smi)
            except:
                mol = None

        if mol:
            smi_list.append(_smi)
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")
    # with open('crossover_nothing.smi', 'w') as f:
    #     for smi in smi_list:
    #         f.write(f'{smi}\n')

    start_time = time.time()
    print("-- start -- our crossover with swap")
    n_valid = 0
    n_fail = 0

    smi_list = []
    for _ in range(100000):
        smi1, smi2 = random.sample(smiles, 2)
        gate = 0
        mol = None
        new_smi = None

        try:
            new_smi, mol = crossover_smiles(smi1, smi2,
                                            prepare_rigid_crossover, True)
        except PermissionError:
            pass
        if not mol:
            try:
                new_smi, mol = crossover_smiles(smi2, smi1,
                                                prepare_rigid_crossover, True)
            except PermissionError:
                pass

        if mol:
            smi_list.append(new_smi)
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")
    # with open('crossover_our_swap.smi', 'w') as f:
    #     for smi in smi_list:
    #         f.write(f'{smi}\n')

    start_time = time.time()
    print("-- start -- our crossover")
    n_valid = 0
    n_fail = 0

    smi_list = []
    for _ in range(100000):
        smi1, smi2 = random.sample(smiles, 2)
        gate = 0
        mol = None
        new_smi = None

        try:
            new_smi, mol = crossover_smiles(smi1, smi2,
                                            prepare_rigid_crossover, True)
        except PermissionError:
            pass

        if mol:
            smi_list.append(new_smi)
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")
    # with open('crossover_our.smi', 'w') as f:
    #     for smi in smi_list:
    #         f.write(f'{smi}\n')

    start_time = time.time()
    print("-- start -- replace atom !")
    n_valid = 0
    n_fail = 0

    for _ in range(100000):
        smi1 = random.sample(smiles, 1)[0]
        try:
            _smi, _mol = replace_atom(smi1)
            mol = Chem.MolFromSmiles(_smi)
        except:
            mol = None

        if mol:
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")

    start_time = time.time()
    print("-- start -- delete atom !")
    n_valid = 0
    n_fail = 0

    for _ in range(100000):
        smi1 = random.sample(smiles, 1)[0]
        try:
            _smi, _mol = delete_atom(smi1)
            mol = Chem.MolFromSmiles(_smi)
        except:
            mol = None

        if mol:
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")

    start_time = time.time()
    print("-- start -- add atom !")
    n_valid = 0
    n_fail = 0

    for _ in range(100000):
        smi1 = random.sample(smiles, 1)[0]
        try:
            _smi, _mol = add_atom(smi1)
            mol = Chem.MolFromSmiles(_smi)
        except:
            mol = None

        if mol:
            n_valid += 1
        else:
            n_fail += 1

        # if _ % 100000 == 0:
        #     print(f"Success rate: {n_valid/(n_valid + n_fail)*100:.3f}% @{_} Time: {(time.time() - start_time)/60:.3f} min")
    print(f"Final Success rate: {n_valid/(n_valid + n_fail)*100:.3f}%")
    print(f"Cost time: {(time.time() - start_time)/60:.3f} min")
