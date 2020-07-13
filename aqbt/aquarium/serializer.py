"""serializer.py.

Methods for serializing Aquarium models into a standardized JSON schema
format
"""
from typing import Dict
from typing import List
from typing import Union

from pydent.base import ModelBase
from pydent.models import FieldValue
from pydent.models import Sample

DataType = Dict[str, Union[str, int, dict]]


def _url_join(*x):
    sep = "/"
    return sep.join([str(_x).strip(sep) for _x in x])


def uri(model: ModelBase):
    session = model.session
    return _url_join(session.url, model.get_tableized_name(), model._primary_key)


def model_serializer(
    model,
    only: Union[str, tuple, dict, list] = None,
    include: Union[str, tuple, dict, list] = None,
) -> DataType:
    if model is None:
        return None
    return model.dump(
        only=only, include=include, include_model_type=True, include_uri=True
    )


def field_value_serializer(fv: FieldValue) -> DataType:
    if fv is None:
        return None
    if fv.field_type.ftype == "sample":
        return {
            "__ftype__": "sample",
            "value": model_serializer(fv.sample, only=("id", "name")),
        }
    else:
        return {"__ftype__": fv.field_type.ftype, "value": fv.value}


def _field_value_properties_serializer(sample) -> dict:
    ft_dict = {}
    for ft in sample.sample_type.field_types:
        ft_dict[ft.name] = ft
    fv_dict = sample.field_value_dictionary()

    properties = {}
    for ftname, ft in ft_dict.items():
        if ft.array:
            properties[ftname] = {
                "__ftype__": "array[{}]".format(ft.ftype),
                "value": [field_value_serializer(fv) for fv in fv_dict[ftname]],
            }
        else:
            properties[ftname] = field_value_serializer(fv_dict[ftname])
    return properties


def sample_serializer(
    sample: Sample,
    only: Union[dict, list, tuple, str] = None,
    include: Union[dict, list, tuple, str] = None,
) -> DataType:
    data = model_serializer(sample, only=only, include=include)
    data["properties"] = _field_value_properties_serializer(sample)
    data["sample_type"] = model_serializer(sample.sample_type, only=("id", "name"))
    return data


# TODO: rename
# TODO: sample_type
def samples_serializer(
    samples: List[Sample],
    only: Union[dict, list, tuple, str] = None,
    include: Union[dict, list, tuple, str] = None,
) -> List[DataType]:
    session = samples[0].session
    data = []
    with session.with_cache(True, using_models=True) as sess:
        g = sess.browser.sample_network(samples)
        for model_name, sample_id in g.nodes:
            if model_name == "Sample":
                data.append(
                    sample_serializer(
                        sess.Sample.find(sample_id), only=only, include=include
                    )
                )
    return data
