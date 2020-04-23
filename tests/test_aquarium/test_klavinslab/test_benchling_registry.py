from aqbt.aquarium import Linter


def test_registry_connector_find(registry, aquarium):
    sample = aquarium.Sample.find(26889)
    seq = registry.connector.find(26889)
    assert seq


def test_registry(registry, aquarium):
    sample = aquarium.Sample.find(26889)
    seq = registry.get_sequence(sample)
    assert seq


def test_fake_registry(registry):
    registry.use_fake_cache(100)


def test_make_pcr_fragment(registry):
    production = registry.session
    fragment = production.Sample.find(32717)
    linter = Linter()
    linter.lint_fragment(registry, fragment)
    print(linter.report())


def test_find(benchling):
    assert benchling.DNASequence.find("seq_Kzxlbux9")


def test_find_by_share_link(benchling):
    assert benchling.DNASequence.from_share_link(
        "https://benchling.com/s/seq-wyfPJ9kHBicShLrsQfzX"
    )


# def test_registry_crawler():
#     session = sessions.aqproduction
#     crawler = RegistryCrawler(session, sessions.klregistry)
#
#     for samples in pydent_utils.sample_pages(
#         session=session,
#         query={"sample_type_id": [d.id for d in crawler.dna_types]},
#         num_pages=5,
#     ):
#         crawler.crawl(samples)


# def test_update_benchling_registry_from_aq():
#     from lobio.klavinslab.pydent_utils import is_keystone_op
#     from pydent.base import ModelBase
#     from pydent.models import Sample
#
#     # TODO: auto_associate data for where it was derived from
#     # TODO: specify benchling registry id
#     # TODO: auto register DNA
#
#     def sample_network(session, samples, reverse=False, g=None):
#         """Build a DAG of :class:`Samples <pydent.models.Sample>` from their.
#
#         :class:`FieldValues <pydent.models.FieldValue>`.
#
#         .. versionadded:: 0.1.5a7
#             method added
#
#         :param samples: list of samples
#         :param reverse: whether to reverse the edges of the final graph
#         :param g: the graph
#         :return:
#         """
#
#         dna_types = session.SampleType.where({'name': ['Plasmid', 'Fragment']})
#         valid_dna_sample_type_ids = [d.id for d in dna_types]
#         OP_SAMPLES_KEY = 'op_samples'
#         OP_KEY = 'operation'
#
#         def dna_like(sample):
#             return isinstance(sample, Sample) and \
#                    sample.sample_type_id in valid_dna_sample_type_ids
#
#         def only_dna_like(models):
#             return [m for m in models if dna_like(m)]
#
#         def get_models(model):
#             for fv in model.field_values:
#                 if fv.sample and issubclass(type(fv.sample), ModelBase):
#                     yield fv.sample, {"field_value_name": fv.name}
#             if hasattr(model, OP_SAMPLES_KEY):
#                 op_samples = getattr(model, OP_SAMPLES_KEY)
#                 for op_sample in op_samples:
#                     print("OKOSDKJFLJSDF")
#                     yield op_sample, {'from_operation': getattr(model, OP_KEY)}
#
#         def op_selector(ops):
#             ops = [op for op in ops if op.status != 'planning']
#             if len(ops) == 1:
#                 return ops[0]
#             elif not ops:
#                 return None
#             else:
#                 return None
#                 # raise ValueError("Ambiguous with {} ops".format(len(ops)))
#
#         def cache_plasmids(session, plasmids):
#             fvs = session.FieldValue.where({
#                 'parent_class': 'Operation',
#                 'child_sample_id': [p.id for p in plasmids],
#                 'role': 'output'
#             })
#             ops = session.browser.get(fvs, {
#                 'operation': {
#                     'operation_type': {},
#                     'field_values': 'sample'
#                 }
#             })['operation']
#
#             grouped = {}
#             for op in ops:
#                 if is_keystone_op(op):
#                     fvs = [fv for fv in op.outputs if fv.child_sample_id in [p.id for p in plasmids]]
#                     for fv in fvs:
#                         grouped.setdefault(fv.child_sample_id, list())
#                         grouped[fv.child_sample_id].append(op)
#
#             selected_grouped = {}
#             for sample_id, ops in grouped.items():
#                 selected_grouped[sample_id] = op_selector(ops)
#
#             for sample_id, op in selected_grouped.items():
#                 if op:
#                     op_samples = only_dna_like([fv.sample for fv in op.inputs])
#                     sample = session.Sample.find(sample_id)
#                     setattr(sample, OP_SAMPLES_KEY, op_samples)
#                     setattr(sample, OP_KEY, op)
#
#         def cache_func(models):
#
#             plasmids = [m for m in models if m.sample_type.name == 'Plasmid']
#             cache_plasmids(session, plasmids)
#             session.browser.get(models, {
#                 "field_values": "sample",
#                 'sample_type': {}}
#             )
#
#         return session.browser.relationship_network(
#             samples, get_models=get_models, cache_func=cache_func, reverse=reverse, g=g
#         )
#
#     import networkx as nx
#
#     registry = sessions.klregistry.copy()
#     registry.use_cache = True
#
#     with sessions.aqproduction.with_cache(timeout=60) as session:
#         dna_types = session.SampleType.where({'name': ['Plasmid', 'Fragment']})
#
#         samples = session.Sample.last(50, query={'sample_type_id': [d.id for d in dna_types], 'user_id': 66})
#         g = sample_network(session, samples)
#
#         # setup graph
#         for n, ndata in g.nodes(data=True):
#             ndata['errors'] = []
#             ndata['sequence'] = None
#             ndata['is_registered'] = False
#
#         # check registration
#         for n in nx.topological_sort(g):
#             ndata = g.nodes[n]
#             preds = g.predecessors(n)
#             sample = session.Sample.find(n[1])
#             assert sample
#
#             is_registered = False
#             found_seq = registry.find_in_cache(sample)
#             if found_seq:
#                 ndata['is_registered'] = True
#                 ndata['sequence'] = found_seq
#             ndata['is_registered'] = is_registered
#
#         # first get the sequence
#         found_using_link = 0
#         for n, ndata in g.nodes(data=True):
#             sample = session.Sample.find(n[1])
#             if not ndata['is_registered']:
#                 seq = registry.get_sequence(sample)
#                 if seq:
#                     ndata['sequence'] = seq
#                     found_using_link += 1
#
#         # check for errors
#         # percolate errors up if sequence is missing
#         # if there are no errors, it is assumed the sequence can be derived
#         for n in nx.topological_sort(g):
#             sample = session.Sample.find(n[1])
#             if sample.sample_type.name == 'Primer':
#                 continue
#             elif sample.sample_type.name not in ['Plasmid', 'Fragment']:
#                 continue
#             ndata = g.nodes[n]
#             if not ndata['sequence']:
#                 parents = list(g.predecessors(n))
#
#                 if not parents:
#                     ndata['errors'].append(
#                         "{} has no predecessors".format(n)
#                     )
#                 else:
#                     parent_errors = []
#                     for parent in parents:
#                         pdata = g.nodes[parent]
#                         parent_errors += pdata['errors']
#                     if parent_errors:
#                         ndata['errors'].append(
#                             "{} parents have errors".format(n)
#                         )
#                     elif sample.sample_type.name == 'Plasmid':
#                         pass
#                     elif sample.sample_type.name == 'Fragment':
#                         if not len(parents) >= 3:
#                             ndata['errors'].append(
#                                 'fragment is missing parents'
#                             )
#
#         # build sequences
#         # set 'sequence'
#         # pull any errors
#         for n in nx.topological_sort(g):
#             ndata = g.nodes[n]
#             print(ndata)
#             sample = session.Sample.find(n[1])
#             if not ndata['sequence'] and not ndata['errors']:
#                 print('building sequence for {}'.format(n))
#                 parents = g.predecessors(n)
#                 if sample.sample_type.name == 'Plasmid':
#                     sequences = [g.nodes[p]['sequence'] for p in parents]
#                     print(sequences)
#                 elif sample.sample_type.name == 'Fragment':
#                     linter = Linter()
#                     products = linter.lint_fragment(registry, sample)
#                     linter.report()
#                     if linter.errors:
#                         ndata['errors'] += linter.errors
#                     elif len(products) == 1:
#                         product = products[0]
#                         ndata['sequence'] = registry.connector.convert(product, to='DNASequence')
#                     elif len(products) > 1:
#                         ndata['errors'].append('More than one product')
#                 print("built")
#                 print(ndata)
#
#         # register sequences
#         for n, ndata in g.nodes(data=True):
#             sample = session.Sample.find(n[1])
#             if not ndata['is_registered'] and ndata['sequence']:
#                 try:
#                     seq = registry.register(sample, ndata['sequence'], overwrite=False,
#                                             do_raise=False)
#                 except Exception as e:
#                     ndata['errors'].append(str(e))
#                 if seq:
#                     ndata['sequence'] = seq
