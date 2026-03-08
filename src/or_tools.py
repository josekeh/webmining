import random
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import cdist
from ortools.linear_solver import pywraplp


def points_sub(points, i):
    """
    De una lista de numeros sustrae el numero i en una copia nueva.
    Se considera que poinla lista no tiene elementos duplicados (por ser lista de indices) y que i se encuentra en la lista.

    Args:
      points: lista de puntos.
      i: elemento a sustraer de la lista

    Returns:
      new: lista nueva similar a points pero sin el elemento i.
    """
    new = points.copy()
    new.remove(i)
    return new

def load_constraints(modelo, list_index_puntos, x_ij, case = "tsp"):
    """
    Carga las restricciones al modelo.

    Args:
      modelo: modelo de ortools al que se le van a cargar las restricciones.
      list_index_puntos: lista de indices de los puntos.
      x_ij: diccionario de variables binarias del modelo.
      case: string que indica el caso para el que se van a cargar las restricciones. Puede ser "tsp" o "vrp".

    Returns:
      modelo: modelo con las restricciones cargadas.
    """
    if case == "tsp":

        # uso n como n-1 ya que arranco de indice 0 en la lista. Por esto en mtz no resto 1 ni en u_i
        n = max(list_index_puntos)
        
        u_i = {i: modelo.IntVar(1,n,'u_'+str(i)) for i in points_sub(list_index_puntos, 0)}
        
        #1
        for j in list_index_puntos:
            modelo.Add(sum(x_ij[i][j] for i in points_sub(list_index_puntos, j)) == 1)

        #2
        for i in list_index_puntos:
            modelo.Add(sum(x_ij[i][j] for j in points_sub(list_index_puntos, i)) == 1)

        #3
        for i in points_sub(list_index_puntos, 0):
            for j in points_sub(points_sub(list_index_puntos, 0), i):
                modelo.Add(u_i[i] - u_i[j] + 1 <= n * (1- x_ij[i][j]))

    elif case == "vrp":
        # Restricciones para VRP
        pass

    return modelo


def calcular_ruta_optima(puntos, costos):

    # Definimos el modelo
    modelo = pywraplp.Solver.CreateSolver('SAT')

    list_index_puntos = list(puntos.index)

    # Definimos variables

    x_ij={i:{j: modelo.IntVar(0,1,'x_'+str(i)+'_'+str(j)) for j in points_sub(list_index_puntos, i)} for i in list_index_puntos}
    # u_i = {i: modelo.IntVar(1,n,'u_'+str(i)) for i in points_sub(list_index_puntos, 0)}

    # Definimos la función objetivo
    obj_expr = sum(x_ij[i][j] * costos[i][j] for i in list_index_puntos for j in points_sub(list_index_puntos, i))

    #Restricciones:

    modelo =  load_constraints(modelo, list_index_puntos, x_ij, case = "tsp")

    # Solver
    modelo.Minimize(obj_expr)

    # Ejecutamos el solver
    status = modelo.Solve()

    return status, modelo, x_ij, list_index_puntos, costos

