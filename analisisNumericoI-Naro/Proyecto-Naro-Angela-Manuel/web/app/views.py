# -*- coding: utf-8
"""contiene las vistas, que son funciones que importan código html
y que están ligadas mediante la función de enrutamiento route()"""

from flask import render_template, flash, redirect
from app import app, numericoUtils, factolu, Cholesky_mod as ch, forms
import numpy
from time import time




A = numpy.array([[1, 1, 1, 1, 1, 1],
                [2, 1, 3, 2, 1, 1],
                [1, 1, 1, 0, 0, 0],
                [1, 1, 0, 0, 1, 1]])

A_str = str(A)

B = A[:, :4]
B_str = str(B)

C = A[:, 4:]
C_str = str(C)

b = numpy.array([[170],
                [320],
                [70],
                [70]])
b_str = str(b)

# Preparando sistema 1
c1 = -numpy.matmul(C, numpy.array([[1], [0]]))
sist_1 = numpy.c_[B, c1]
sist_1_str = str(sist_1)

# Preparando sistema 2
c2 = -numpy.matmul(C, numpy.array([[0], [1]]))
sist_2 = numpy.c_[B, c2]
sist_2_str = str(sist_2)

# Preparando sistema 3
sist_3 = numpy.c_[B, b]
sist_3_str = str(sist_3)

opcion = -1


@app.route('/')
@app.route('/index')
def index():
    """renderiza la página principal"""
    return render_template("index.html")

@app.route('/inicio')
def inicio():
    """renderiza la página inicio.html"""
    return render_template("inicio.html")

@app.route('/problema')
def problema():
    """renderiza la página que presenta el problema"""
    return render_template("problema.html")

@app.route('/planteamiento', methods = ['GET', 'POST'])
def planteamiento():
    posts = {'A': A_str,
            'b': b_str,
            'B': B_str,
            'C': C_str,
            'sistema1': sist_1_str,
            'sistema2': sist_2_str,
            'sistema3': sist_3_str}
    form = forms.method_form()
    if form.validate_on_submit() :
        flash(form.metodo.data)
        global opcion
        opcion = form.metodo.data
        return redirect('/solucion')
    return render_template("planteamiento.html", posts = posts, form = form)

@app.route('/solucion')
def solucion():
    """renderiza la página que presenta la solución del problema"""
    global opcion
    if opcion == 'elim_Gauss' :
        # Sistema 1
        t_inicio = time()
        x1 = numericoUtils.elimGauss(B, c1);
        t_final = time()
        x1_str = str(x1)
        t_p = t_final - t_inicio
        t1_str = str(t_p)
        # Sistema 2
        t_inicio = time()
        x2 = numericoUtils.elimGauss(B, c2);
        t_final = time()
        x2_str = str(x2)
        t_p = t_final - t_inicio
        t2_str = str(t_p)
        # Sistema 3
        t_inicio = time()
        x3 = numericoUtils.elimGauss(B, b);
        t_final = time()
        x3_str = str(x3)
        t_p = t_final - t_inicio   
        t3_str = str(t_p)
    elif opcion == 'factolu' :
        t_inicio = time()
        x1 = factolu.sol(B, c1);
        t_final = time()
        x1_str = str(x1)
        t_p = t_final - t_inicio
        t1_str = str(t_p)
        # Sistema 2
        t_inicio = time()
        x2 = factolu.sol(B, c2);
        t_final = time()
        x2_str = str(x2)
        t_p = t_final - t_inicio
        t2_str = str(t_p)
        # Sistema 3
        t_inicio = time()
        x3 = factolu.sol(B, b);
        t_final = time()
        x3_str = str(x3)
        t_p = t_final - t_inicio   
        t3_str = str(t_p)
    elif opcion == 'cholesky' :
        if ch.posdef(B) == 0 :
            t_inicio = time()
            x1 = ch.choleskyporTrans(numpy.mat(B), numpy.mat(c1))
            t_final = time()
            x1_str = str(x1)
            t_p = t_final - t_inicio
            t1_str = str(t_p)

            t_inicio = time()
            x2 = ch.choleskyporTrans(numpy.mat(B), numpy.mat(c2))
            t_final = time()
            x2_str = str(x2)
            t_p = t_final - t_inicio
            t2_str = str(t_p)

            t_inicio = time()
            x3 = ch.choleskyporTrans(numpy.mat(B), numpy.mat(b))
            t_final = time()
            x3_str = str(x3)
            t_p = t_final - t_inicio
            t3_str = str(t_p)
        else:
            if ch.simetric(B) == 0:
                t_inicio = time()
                x1 = ch.choleskyporTrans(B, c1)
                t_final = time()
                x1_str = str(x1)
                t_p = t_final - t_inicio
                t1_str = str(t_p)

                t_inicio = time()
                x2 = ch.choleskyporTrans(B, c2)
                t_final = time()
                x2_str = str(x)
                t_p = t_final - t_inicio
                t2_str = str(t_p)

                t_inicio = time()
                x3 = ch.choleskyporTrans(B, b)
                t_final = time()
                x3_str = str(x)
                t_p = t_final - t_inicio
                t3_str = str(t_p)
            else:
                t_inicio = time()
                x1 = ch.choleskyTrans(B, c1)
                t_final = time()
                x1_str = str(x1)
                t_p = t_final - t_inicio
                t1_str = str(t_p)

                t_inicio = time()
                x2 = ch.choleskyTrans(B, c2)
                t_final = time()
                x2_str = str(x2)
                t_p = t_final - t_inicio
                t2_str = str(t_p)

                t_inicio = time()
                x3 = ch.choleskyTrans(B, b)
                t_final = time()
                x3_str = str(x3)
                t_p = t_final - t_inicio
                t3_str = str(t_p)

    posts = [ 
            {'sistema': {'x': x1_str, 't': t1_str}, 'nombre': 'sistema 1'},
            {'sistema': {'x': x2_str, 't': t2_str}, 'nombre': 'sistema 2'},
            {'sistema': {'x': x3_str, 't': t3_str}, 'nombre': 'sistema 3'}
            ]
    return render_template("solucion.html", posts = posts)


@app.route('/creditos')
def creditos():
    """renderiza la página de creditos"""
    return render_template("creditos.html")
