import networkx as nx
from typing import Hashable, Generator, TypeVar, List, Tuple, Optional, Dict, Callable
from pydent import AqSession
from pydent import ModelBase
from tqdm.auto import tqdm
import pandas as pd

Node = Hashable
_T = TypeVar('_T')


def is_leaf(n: Node, g: nx.DiGraph):
    successors = list(g.successors(n))
    if successors:
        return False
    return True


def is_root(n: Node, g: nx.DiGraph):
    predecessors = list(g.predecessors(n))
    if predecessors:
        return False
    return True


def iter_leaves(g: nx.DiGraph) -> Generator[Node, None, None]:
    for n in g.nodes():
        if is_leaf(n, g):
            yield n


def iter_roots(g: nx.DiGraph) -> Generator[Node, None, None]:
    for n in g.nodes():
        if is_root(n, g):
            yield n


def iter_split_sets(arr: List[_T]) -> Tuple[List[_T], _T]:
    for i in range(len(arr)):
        yield tuple(arr[:i] + arr[i + 1:]), arr[i]


class HashFunctions(object):

    @staticmethod
    def int_list_hash(ids: List[int], prefix: Optional[str] = None) -> Tuple[str, Tuple[int]]:
        return prefix, tuple(sorted(ids))


def group_by_sample_id(items):
    data = {}
    for item in items:
        data.setdefault(item.sample_id, list())
        data[item.sample_id].append(item)
    return data


def group_by_object_type_id(items):
    data = {}
    for item in items:
        data.setdefault(item.object_type_id, list())
        data[item.object_type_id].append(item)
    return data


def is_aq_model(x):
    return issubclass(x.__class__, ModelBase)


def get_models(sess, x, model):
    interface = sess.model_interface(model)
    if isinstance(x, int):
        res = [interface.find(x)]
        if res[0] is None:
            raise ValueError("Could not find {}".format(x))
    elif isinstance(x, str):
        res = [interface.find_by_name(x)]
        if res[0] is None:
            raise ValueError("Could not find {}".format(x))
    elif is_aq_model(x):
        res = [x]
        if res is None:
            raise ValueError
    elif isinstance(x, dict):
        res = interface.where(x)
        if res is None:
            raise ValueError
    elif isinstance(x, (list, tuple)):
        res = []
        for _x in x:
            res.extend(get_models(sess, _x, model))
    else:
        raise TypeError
    return res


class LIMSDB(object):
    YEAST_OBJECT_TYPES = ['Yeast Glycerol Stock', 'Yeast Plate', 'Yeast Competent Aliquot', 'Yeast Competent Cell']
    DNA_OBJECT_TYPES = ['Plasmid Stock', 'Plasmid Glycerol Stock', '1 ng/µL Plasmid Stock', 'Fragment Stock',
                        '1 ng/µL Fragment Stock']

    def __init__(self, sess: AqSession, strains=None, n=None, timeout: int = 60, progressbar: Callable = tqdm):
        """

        :param sess:
        :param strains:
        :param n:
        :param timeout:
        :param progressbar:
        """
        self.progressbar = progressbar
        self.sess = sess.with_cache(timeout=timeout)
        self.common = {
            'sample_types': {
                'yeast': self.sess.SampleType.find_by_name('Yeast Strain'),
                'plasmid': self.sess.SampleType.find_by_name('Plasmid'),
                'fragment': self.sess.SampleType.find_by_name('Fragment'),
            },
            'object_types': {
                'yeast': get_models(self.sess, self.YEAST_OBJECT_TYPES, 'ObjectType'),
                'dna': get_models(self.sess, self.DNA_OBJECT_TYPES, 'ObjectType')
            }
        }
        self.prov: nx.DiGraph = self.build_prov_graph(strains=strains, n=n)
        self.cache_items()

    def node_to_hash_fn(self, x, prefix=None):
        node, ndata = x
        return self.hash_fn(ndata['integrants'], prefix=prefix)

    def find_by_hash(self, x, prefix=None) -> List[Tuple[Node, Dict]]:
        k1 = self.hash_fn(x, prefix=prefix)
        visited = []
        for x in self.prov.nodes(data=True):

            ndata = x[1]
            if 'mating_type' in ndata:
                k2 = self.node_to_hash_fn(x, prefix=ndata['mating_type'])
                if k1 == k2:
                    visited.append(x)
        return visited

    def hash_fn(self, x, prefix=None):
        return HashFunctions.int_list_hash(x, prefix=prefix)

    def build_prov_graph(self, strains: Optional[List[ModelBase]] = None, n=None) -> nx.DiGraph:
        g = self._build_prov_graph(strains, n=n)
        grev = g.reverse()
        for n in grev.nodes():
            integrants = list(self._get_integrants(n, grev))
            g.nodes[n]['integrants'] = tuple(sorted(integrants))
        return g

    @staticmethod
    def _get_integrants(n, grev):
        for n1, n2 in list(nx.dfs_edges(grev, n)):
            edata = grev[n1][n2]
            if edata['field_name'] == 'Integrant':
                yield n2

    def _build_prov_graph(self, strains: Optional[List[ModelBase]] = None, n=None) -> nx.DiGraph:
        g = nx.DiGraph()

        def add_sample_node(s, **kwargs):
            g.add_node(s.id, sample_type=s.sample_type_id, sample=s.dump(), **kwargs)

        if not strains:
            assert n
            strains = self.sess.Sample.last(n, query={
                'sample_type_id': self.common['sample_types']['yeast'].id,
                'user_id': 66
            })
        while strains:
            if strains:
                self.sess.browser.get(strains, {
                    'field_values': {
                        'field_type': 'allowable_field_types',
                        'sample': []
                    },
                    'sample_type': {
                        'field_types': 'allowable_field_types'
                    }
                })

                parents = {}
                for s in strains:
                    add_sample_node(s, mating_type=s.properties['Mating Type'])

                    parent = None
                    if 'Parent' in s.properties:
                        parent = s.properties['Parent']

                    integrant = None
                    if 'Integrant' in s.properties:
                        integrant = s.properties['Integrant']

                    if parent is not None:
                        parents[parent.id] = parent
                        add_sample_node(parent, mating_type=parent.properties['Mating Type'])
                        g.add_edge(parent.id, s.id, field_name='Parent')

                    if integrant is not None:
                        add_sample_node(integrant)
                        g.add_edge(integrant.id, s.id, field_name='Integrant', sample=integrant.dump())

            new_parents = {}
            sids = [s.id for s in strains]
            for sid, p in parents.items():
                if sid not in sids:
                    new_parents[sid] = p
            strains = list(new_parents.values())
        return g

    def find_in_cache(self, query: dict, model_name: str):
        visited = []
        for item in self.sess.browser.model_cache[model_name].values():
            for k, v in query.items():
                if getattr(item, k) == v:
                    visited.append(item)
        return visited

    def find_yeast_items(self, sample_id):
        items = self.find_items(sample_id, [ot.id for ot in self.common['object_types']['yeast']])
        return items

    def find_dna_items(self, sample_id):
        items = self.find_items(sample_id, [ot.id for ot in self.common['object_types']['dna']])
        return items

    def find_items(self, sample_id, otids):
        visited = []
        for item in self.sess.browser.model_cache["Item"].values():
            if item.sample_id == sample_id and item.location != 'deleted' and item.object_type_id in otids:
                visited.append(item)
        return visited

    def cache_items(self):
        all_ots = []
        for k, v in self.common['object_types'].items():
            all_ots.extend(v)
        ot_by_st = {}
        for ot in all_ots:
            ot_by_st.setdefault(ot.sample_type_id, list())
            ot_by_st[ot.sample_type_id].append(ot)

        for stid, ots in ot_by_st.items():
            samples = self.find_in_cache({'sample_type_id': stid}, 'Sample')
            self._cache_items([s.id for s in samples], ots)

    def _cache_items(self, sample_ids: List[int], object_types):
        sess = self.sess
        ots = get_models(sess, object_types, 'ObjectType')
        items = []
        for ot in self.progressbar(ots):
            query = {'sample_id': sample_ids, 'object_type_id': ot.id}
            items_ = sess.Item.where(query)
            items_ = [i for i in items_ if i.location != 'deleted']
            items.extend(items_)
        return items

    def create_trajectories(self, prefix, parts_arr, g=None) -> nx.DiGraph:
        """
        Create all possible build trajectories.

        :param prefix:
        :param parts_arr:
        :param g:
        :return:
        """
        if g is None:
            g = nx.DiGraph()
        for parts_arr_ in parts_arr:
            n1 = self.hash_fn(prefix, parts_arr_)
            visited = []
            for parts, x in iter_split_sets(parts_arr_):
                if len(parts) > 1:
                    visited.append(parts)
                n2 = self.hash_fn(prefix, parts)
                g.add_edge(n2, n1, plasmid_id=x)
            self.create_trajectories(prefix, visited, g=g)
        return g


def by_key_value(arr, keyfn, valuefn, iffn=None):
    data = {}
    for a in arr:
        key = keyfn(a)
        data.setdefault(key, list())
        v = valuefn(a)
        if iffn and iffn(v):
            data[key].append(valuefn(a))
    return data


from copy import deepcopy


class BuildPaths(object):

    def __init__(self, db: LIMSDB):
        self.db = db

    def create_trajectories(self, prefix, parts_arr, g=None) -> nx.DiGraph:
        """
        Create all possible build trajectories.

        :param prefix:
        :param parts_arr:
        :param g:
        :return:
        """
        hash_fn = self.db.hash_fn
        if g is None:
            g = nx.DiGraph()
        for parts_arr_ in parts_arr:
            n1 = hash_fn(x=parts_arr_, prefix=prefix)
            visited = []
            for parts, x in iter_split_sets(parts_arr_):
                if len(parts) > 1:
                    visited.append(parts)
                n2 = hash_fn(x=parts, prefix=prefix)
                g.add_edge(n2, n1, plasmid_id=x)
            self.create_trajectories(prefix, visited, g=g)
        return g

    @staticmethod
    def validate_parts(parts):

        for part in parts:
            assert 'gate' in part
            assert 'haploid' in part
            assert 'sample' in part

    def generate_build_paths(self, parts: List[dict]):
        parts = deepcopy(parts)
        for part in parts:
            sample = self.db.sess.Sample.find_by_name(part['sample']['name'])
            part['sample'] = sample.dump()
        self.validate_parts(parts)
        partsdict = by_key_value(parts, lambda x: x['haploid'], lambda x: x['sample']['id'], lambda x: True)
        trajectories = {}
        print(partsdict)
        for haploid, parts in partsdict.items():
            traj = self.create_trajectories(prefix=haploid, parts_arr=[parts])

            nodelist = list(traj.nodes())
            traj.add_node("START")

            for n1 in nodelist:
                prefix, parts_ = n1
                if prefix == 'Mat A':
                    prefix = 'MATa'
                elif prefix == 'Mat Alpha':
                    prefix = 'MATalpha'
                else:
                    raise RuntimeError('prefix "{}" unexpected for {}'.format(prefix, n1))

                found = self.db.find_by_hash(parts_, prefix=prefix)
                for sample_id, ndata2 in found:
                    items = self.db.find_yeast_items(sample_id)
                    if items:
                        traj.add_edge('START', n1, weight=0)

            for n1, n2, edata in traj.edges(data=True):
                if n1 != 'START':
                    plasmid_id = edata['plasmid_id']
                    items = self.db.find_dna_items(plasmid_id)
                    if items:
                        edata['weight'] = 1
                    else:
                        edata['weight'] = 4

            trajectories[haploid] = traj
        return trajectories


def by_key_value(arr, keyfn, valuefn, iffn=None):
    data = {}
    for a in arr:
        key = keyfn(a)
        data.setdefault(key, list())
        v = valuefn(a)
        if iffn and iffn(v):
            data[key].append(valuefn(a))
    return data

class BuildPaths(object):

    def __init__(self, db: LIMSDB):
        self.db = db

    def create_trajectories(self, prefix, parts_arr, g=None) -> nx.DiGraph:
        """
        Create all possible build trajectories.

        :param prefix:
        :param parts_arr:
        :param g:
        :return:
        """
        hash_fn = self.db.hash_fn
        if g is None:
            g = nx.DiGraph()
        for parts_arr_ in parts_arr:
            n1 = hash_fn(x=parts_arr_, prefix=prefix)
            visited = []
            for parts, x in iter_split_sets(parts_arr_):
                if len(parts) > 1:
                    visited.append(parts)
                n2 = hash_fn(x=parts, prefix=prefix)
                g.add_edge(n2, n1, plasmid_id=x)
            self.create_trajectories(prefix, visited, g=g)
        return g

    @staticmethod
    def validate_parts(parts):

        for part in parts:
            assert 'gate' in part
            assert 'haploid' in part
            assert 'sample' in part



    def generate_build_paths(self, parts: List[dict]):
        parts = deepcopy(parts)
        for part in parts:
            sample = self.db.sess.Sample.find_by_name(part['sample']['name'])
            part['sample'] = sample.dump()
        self.validate_parts(parts)
        partsdict = by_key_value(parts, lambda x: x['haploid'], lambda x: x['sample']['id'], lambda x: True)
        trajectories = {}
        print(partsdict)
        for haploid, parts in partsdict.items():
            traj = self.create_trajectories(prefix=haploid, parts_arr=[parts])

            nodelist = list(traj.nodes())
            traj.add_node("START")

            for n1 in nodelist:
                prefix, parts_ = n1
                if prefix == 'Mat A':
                    prefix = 'MATa'
                elif prefix == 'Mat Alpha':
                    prefix = 'MATalpha'
                else:
                    raise RuntimeError('prefix "{}" unexpected for {}'.format(prefix, n1))

                found = self.db.find_by_hash(parts_, prefix=prefix)
                for sample_id, ndata2 in found:
                    items = self.db.find_yeast_items(sample_id)
                    if items:
                        traj.add_edge('START', n1, weight=0)

            for n1, n2, edata in traj.edges(data=True):
                if n1 != 'START':
                    plasmid_id = edata['plasmid_id']
                    items = self.db.find_dna_items(plasmid_id)
                    if items:
                        edata['weight'] = 1
                    else:
                        edata['weight'] = 4

            trajectories[haploid] = traj
        return trajectories

    def generate_instructions(self, parts):
        self.validate_parts(parts)

        instructions = []

        for h, g in self.generate_build_paths(parts).items():
            leaves = list(iter_leaves(g))
            assert len(leaves) == 1
            paths = list(nx.all_shortest_paths(g, 'START', leaves[0]))
            path_weights = []
            for path in paths:
                weights = []
                for n1, n2 in nx.utils.pairwise(path):
                    weights.append(g[n1][n2]['weight'])
                path_weights.append(tuple(weights))
            sorted_paths = sorted(list(zip(paths, path_weights)), key=lambda x: x[1])
            best_path = sorted_paths[0]
            for i, (n1, n2) in enumerate(nx.utils.pairwise(best_path[0])):
                if n1 == "START":
                    continue
                p1 = n1[0]
                p2 = n2[0]
                _prefix = {
                    'Mat A': 'MATa',
                    'Mat Alpha': 'MATalpha'
                }
                p1 = _prefix[p1]
                p2 = _prefix[p2]

                s1 = self.db.find_by_hash(list(n1[1]), prefix=p1)
                s2 = self.db.find_by_hash(list(n2[1]), prefix=p2)
                edata = g[n1][n2]

                if s1:
                    parent = [(_s[0], _s[1]['sample']['name']) for _s in s1]
                else:
                    parent = n1

                if s2:
                    result = [(_s[0], _s[1]['sample']['name']) for _s in s2]
                else:
                    result = n2

                instructions.append({
                    'step': i,
                    'parent': parent,
                    'dna': edata['plasmid_id'],
                    'result': result,
                    'weight': edata['weight']
                })

        return instructions

    def parse_builds(self, builds):
        all_rows = []
        for build in builds:
            parts = build['parts']
            data = dict(build)
            data['id'] = data['name'] + "_" + str(data['permutation'])
            del data['parts']
            rows = self.generate_instructions(parts)
            for row in rows:
                row.update(data)
            all_rows.extend(rows)
        return pd.DataFrame(all_rows)