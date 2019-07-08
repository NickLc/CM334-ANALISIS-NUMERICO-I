import Tkinter
import time
import tkMessageBox
#root.state(newstate="normal") normal
#root.state(newstate="withdraw") oculto
#root.state(newstate="iconic") minimizar
#
#
#
def aparecer():
    root.state(newstate='iconic')
    time.sleep(5)
    root.state(newstate='normal')
def valor():
    if entry.get() == "Naro":
        print "Bienvenido Naro e.e"
    else:
        print "Tu no eres naro, fuera oh ctm"
def panel():
    if usr.get() == "admin" and pwd.get() == "admin":
        tkMessageBox.showinfo("Genial", "Bienvenido admin :3")
    else:
        tkMessageBox.showerror("Error", "Vete oh mrd! XD")

#Creando Ventana#
root = Tkinter.Tk()
root.title("Panel GSINT")
root.geometry("250x150+500+250")

#Label1
etiqueta1=Tkinter.Label(root, text="Usuario")#, width=100, height=100, anchor="center")
etiqueta1.pack()

#Entry1
usr = Tkinter.Entry(root)
usr.pack()

#Label2
etiqueta2=Tkinter.Label(root, text="Password")#, width=100, height=100), anchor="center")
etiqueta2.pack()

#Entry2
pwd = Tkinter.Entry(root, show="*")
pwd.pack()

#Boton
boton=Tkinter.Button(root, text="Ingresar", command=panel)
boton.pack()


root.mainloop()

