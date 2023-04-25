#!/bin/env python3

from typing import Optional

# the singleton class to print any values.
class Printer:
    # the singleton instance.
    __instance: Optional["Printer"] = None


    # the static method to get the singleton instance.
    @staticmethod
    def get_instance() -> Optional["Printer"]:
        if Printer.__instance is None:
            Printer()
        return Printer.__instance


    # the private constructor.
    def __init__(self) -> None:
        if Printer.__instance is not None:
            raise Exception("This class is a singleton!")
        else:
            Printer.__instance = self


    # this method is to print any values.
    def print(self, value: object) -> None:
        print(value)


def main() -> None:
    try:
        # インスタンス化される
        printer = Printer()
        printer.print(10)

        # インスタンス化されず、例外が発生する
        printer2 = Printer()
        printer2.print("Hello, world!")
    except Exception as e:
        print(e)


if __name__ == "__main__":
    main()