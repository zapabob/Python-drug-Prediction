#!/bin/env python3

from typing import Optional
import tkinter as tk
import classes as cl
import pandas as pd
import time
import threading
import random

# The singleton class to display the main window.
class MainWindow:
    # The singleton instance.
    __instance: Optional["MainWindow"] | None = None


    # The static method to get the singleton instance.
    @staticmethod
    def get_instance(title: str) -> Optional["MainWindow"]:
        if MainWindow.__instance is None:
            MainWindow(title)
        return MainWindow.__instance


    # The private constructor.
    # This constructor receives the window title of string.
    def __init__(self, title: str) -> None:
        if MainWindow.__instance is not None:
            raise Exception("This class is a singleton!")
        else:
            MainWindow.__instance = self

        # set the window title
        self.title: str = title

        # make the window
        self.make_window()

        # run the main loop
        self.window.mainloop()


    # this mathod is to make the window.
    def make_window(self) -> None:
        # initialize the window
        self.window = tk.Tk()

        # set the title of the window
        self.window.title(self.title)

        # make the text box and set to the window
        self.text_box = tk.Entry(self.window, width=40)
        self.text_box.pack(pady=10)

        # make the button and set to the window
        self.button = tk.Button(self.window, text='Get ADME Data', command=self.button_click)
        self.button.pack(pady=10)

        # make the text box to display the result and set to the window.
        # This text box is initially disabled.
        self.result_text = tk.Text(self.window, width=80, height=20, state='disabled')
        self.result_text.pack(pady=10)


    # The method to be called when the button is clicked.
    def button_click(self) -> None:
        # disable the button and change the text of the button
        self.button.configure(state='disabled', text='Getting ADME Data...')

        # get the IUPAC name from the text box
        iupac: str = self.text_box.get()
        smiles: str = cl.SubstanceName(iupac=iupac).get()
        adme_data: pd.DataFrame = cl.SwissADME(smiles).get()

        # display results of adme_data in the text box
        self.result_text.configure(state='normal')
        self.result_text.delete(1.0, tk.END)
        text: str = ""
        for _, item in adme_data.iterrows():
            text = str(item)
        self.result_text.insert(tk.END, text)
        self.result_text.configure(state='disabled')

        # launch the thread to activate the button after random seconds in 1 to 30.
        dulation: int = random.randint(1, 30)
        threading.Thread(target=self.activate_button, args=(dulation,)).start()


    # the method of the thread to activate the button
    # arguments:
    # - dulation (int): the time to wait in seconds
    def activate_button(self, dulation: int) -> None:
        # change button text to count down
        for i in range(dulation, 0, -1):
            self.button.configure(text=f"{i}")
            time.sleep(1)
        
        # change the button text to the initial value and change tht button state to normal
        self.button.configure(text="Get ADME Data")
        self.button.configure(state='normal')
