# -*- coding: utf-8 -*-
"""contiene la estructura predefinida de los forms del html"""

import flask_wtf
import wtforms

class method_form(flask_wtf.Form):
    """esquema del formulario para escoger el método de resolución"""
    metodo = wtforms.RadioField('metodo', choices = [
                                                ('elim_Gauss', 'Eliminación de Gauss'),
                                                ('factolu', 'Factorización LU'),
                                                ('cholesky', 'Cholesky')],
                                                default = False,
                                                validators = [wtforms.validators.DataRequired()])
