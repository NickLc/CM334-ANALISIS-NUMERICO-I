# -*- coding: utf-8 -*-
#!flask/bin/python

from app import app

app.run(host = '127.0.0.1', port = 8081, debug = True)