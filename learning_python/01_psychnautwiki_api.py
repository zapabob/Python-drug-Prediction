#!/bin/env python3

import requests
import re


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

    def get_result(self) -> str:
        return self.result


def main() -> None:
    # make the instance of the class
    psychonautwiki_api = PsychonautWikiAPI(substance_name="methamphetamine")

    # get the result
    result: str = psychonautwiki_api.get_result()

    # print the result
    print(result)


if __name__ == "__main__":
    main()