#!/bin/env python3

import requests
import re
import json
from typing import Any
import sys


class PostJson:
    def __init__(self, url: str, json: dict[str, str], headers: dict[str, str] | None = None) -> None:
        self.url = url
        self.json = json
        self.headers = headers
        self.return_value: str = ""

    def post(self) -> None:
        response = requests.post(self.url, json=self.json, headers=self.headers)
        self.return_value = response.text

    def get_value(self) -> str:
        return self.return_value


class PsychonautWikiAPI:
    def __init__(self, substance_name: str) -> None:
        self.substance_name = substance_name
        self.query = self.get_query()
        self.json: dict[str, str] = {"query": self.query}
        self.headers: dict[str, str] = {
            "Accept-Encoding": "gzip, deflate, br",
            "Content-Type": "application/json",
            "Connection": "keep-alive",
            "DNT": "1"
        }

        # POST
        self.postjson = PostJson(url="https://api.psychonautwiki.org/", 
                                 json=self.json, 
                                 headers=self.headers)
        self.postjson.post()

        # get the response
        self.result: str = self.postjson.get_value()

    def get_query(self) -> str:
        ret: str = ""
        with open("01_query.txt") as f:
            for line in f.readlines():
                line = re.sub(r"<substance_name>", self.substance_name, line)
                ret += line
        return ret

    def get_result(self) -> Any:
        jsondata = json.loads(self.result)

        if jsondata["data"]["substances"] == []:
            raise Exception("The substance name is not found.")

        return jsondata


def main() -> None:
    # get arguments
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <substance_name> [json_file_name]")
        return
    else:
        substance_name = sys.argv[1]
        if len(sys.argv) == 3:
            json_file_name = sys.argv[2]
        else:
            json_file_name = ""

    # PsychonautWikiよりデータを取得
    try:
        result: Any = PsychonautWikiAPI(substance_name).get_result()
    except Exception as e:
        print(e)
        return

    # save the json data
    if json_file_name != "":
        with open(json_file_name, "w") as f:
            json.dump(result, f, indent=4)
            print(f"Saved the json data to {json_file_name}.")
    else:
        print(result)


if __name__ == "__main__":
    main()