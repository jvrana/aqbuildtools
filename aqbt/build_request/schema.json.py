{
  "$schema": "http://json-schema.org/draft-07/schema#",
    "title": "Part JSON",
    "description": "",
    "properties": {
        "name": {
            "type": "string"
        },
        "collection": {
            "type": "string"
        },
        "type": {
            "type": "string",
            "enum": ["composite part"]
        },
        "description": {
            "type": "string"
        },
        "parts": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1
        }
    },
    "required": ["name", "collection", "type", "parts", "description"]
}

composite parts must contain parts while basic parts must contain sequence