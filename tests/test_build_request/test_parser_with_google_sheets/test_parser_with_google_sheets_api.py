import pickle
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import join

from google.auth.transport.requests import Request
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build


SCOPES = ["https://www.googleapis.com/auth/spreadsheets.readonly"]
here = abspath(dirname(__file__))


def get_creds():
    creds = None

    token_path = join(here, "token.pickle")
    creds_path = join(here, "credentials.json")

    if exists(token_path):
        with open(token_path, "rb") as token:
            creds = pickle.load(token)

    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(creds_path, SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open(token_path, "wb") as token:
            pickle.dump(creds, token)
    return creds


def test_main():
    creds = get_creds()
    service = build("sheets", "v4", credentials=creds)

    request = (
        service.spreadsheets()
        .values()
        .get(
            spreadsheetId="1aDZv-A8JRGZbtVFugOsB6w-001pni8DDIoj1IXdRcMQ", range="Design"
        )
    )

    data = request.execute()

    from aqbt.build_request import parse_parts, CellValue

    values = CellValue.to_cell_values(data["values"])

    parsed = parse_parts(values)
    for p in parsed:
        print(p)

    # TODO: create a graph, find leaves, find all intermediate strains
