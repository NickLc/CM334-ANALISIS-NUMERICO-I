# -*- coding: utf-8 -*-
"""Contiene la interfaz gráfica del proyecto."""

import tkinter

class App(tkinter.Frame):
    def __init__(self, master = None):
        super().__init__(master)
        self.master.maxsize(400, 600)
        self.pack()
        self.create_widgets()
    
    def create_widgets(self):
        self.hola = tkinter.Button(self, text = "Hola!", command = self.saludar)
        self.hola.pack(side = "top")
        self.salir = tkinter.Button(self, text = "salir", bg = "red", command = self.master.destroy)
        self.salir.pack(side = "bottom")

    def saludar(self):
        print("Si tuviste un día, mañana será peor!")
    

if __name__ == '__main__' :
    main_w = tkinter.Tk()
    main_w.title('Panel principal')
    main_w.geometry("300x400")
    app1 = App(master = main_w)
    app1.mainloop()
    
