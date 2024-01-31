import scipy
import matplotlib.pyplot as plt
from typing import Union

class ChartParametrs(object):
    """Class ChartParametrs is used to find parameters of flow models

    Attributes
    ----------
    deformations --- [list or tuple] containing full deformations (epsilon)
    
    tensions --- [list or tuple] containing tensions (sigma)
    
    info --- [dict] information about diagramm containing keys: 
        Transition point load
        Transition point unload
        End of linear section point load
        End of linear section point unload
        Limit strength point load
        Limit strength point unload
        Elastic modulus load
        Elastic modulus load
    """

    def __init__(self,
                 deformations: Union[list, tuple],
                 tensions: Union[list, tuple]) -> None:
        try:
            if not deformations or not tensions:
                raise ValueError
            elif len(deformations) != len(tensions):
                raise ValueError
        except ValueError:
            print("TypeError: List or tuples must not be empty,"\
                  " also size of deformations and tensions must"\
                  " be equal.")
        else:
            self.deformations: Union[list, tuple] = deformations
            self.tensions: Union[list, tuple] = tensions
            self.info: dict = dict()

    def __str__(self) -> str:
        pass

    def _find_index_transition_point_load(self) -> int:
        """
        """
        max_index = 0

        for i in range(len(self.deformations)):
            if self.deformations[i] >= self.deformations[max_index]:
                max_index = i

        return max_index

    def _find_index_transition_point_unload(self) -> int:
        """
        """
        index_transition_load = self._find_index_transition_point_load()
        min_index = len(self.deformations) - 1

        for i in range(index_transition_load, len(self.deformations)):
            if self.deformations[i] <= self.deformations[min_index]:
                min_index = i
        
        return min_index

    def _find_index_linear_section_load(self) -> int:
        """
        """
        index_transition_load = self._find_index_transition_point_load()
        x, y = self.deformations, self.tensions
        speed_r_square = 0
        best_speed_r_square = float("inf")
        previous_r_square = 0
        best_index = 0

        for i in range(2, index_transition_load)[::-1]:
            # обработать ошибки линейной регрессии
            parametrs = scipy.stats.linregress(x[:i], y[:i])
            r_square = parametrs[2]
            speed_r_square = abs(r_square - previous_r_square)
            previous_r_square = r_square
            if best_speed_r_square > speed_r_square:
                best_index = i
                best_speed_r_square = speed_r_square

        return best_index
    
    def _find_index_linear_section_unload(self) -> int:
        """
        """
        index_transition_load = self._find_index_transition_point_load()
        index_transition_unload = self._find_index_transition_point_unload()
        x = self.deformations
        y = self.tensions
        speed_r_square = 0
        best_speed_r_square = float("inf")
        previous_r_square = 0
        best_index = 0

        for i in range(index_transition_load + 1, index_transition_unload)[::-1]:
            # обработать ошибки линейной регрессии
            parametrs = scipy.stats.linregress(x[index_transition_load:i], y[index_transition_load:i])
            r_square = parametrs[2]
            speed_r_square = abs(r_square - previous_r_square)
            previous_r_square = r_square
            if best_speed_r_square > speed_r_square:
                best_index = i
                best_speed_r_square = speed_r_square

        return best_index
    
    def _find_limit_strength_load(self, 
                             linear_section_index: int,
                             transition_index: int,
                             admission: float=0.02):
        """
        """
        x = self.deformations[linear_section_index:transition_index]
        y = self.tensions[linear_section_index:transition_index]
        k = self.tensions[linear_section_index] / self.deformations[linear_section_index]

        line_section = lambda x: k * (x - admission)
        interpolation = scipy.interpolate.interp1d(x, y, bounds_error=False)
        # Поправить invalid value encountered in scalar divide
        root = float(scipy.optimize.fsolve(lambda x: line_section(x) - interpolation(x), x[0]))

        return root, line_section(root)

    def _find_limit_strength_unload(self, 
                             linear_section_index: int,
                             transition_index_load: int,
                             transition_index_unload: int,
                             admission: float=0.02):
        """
        """
        x = self.deformations[linear_section_index:transition_index_unload]
        y = self.tensions[linear_section_index:transition_index_unload]
        k = (self.tensions[transition_index_load] - self.tensions[linear_section_index]) / (self.deformations[transition_index_load] - self.deformations[linear_section_index])

        linear_section = lambda x: k * (x + admission - self.deformations[linear_section_index]) + self.tensions[linear_section_index]
        interpolation = scipy.interpolate.interp1d(x, y, bounds_error=False)
        # Поправить invalid value encountered in scalar divide
        # не работает, если начальная точка x[0]
        root = float(scipy.optimize.fsolve(lambda x: linear_section(x) - interpolation(x), x[-1]))

        return root, linear_section(root)

    def _return_elastic_modulus_load(self,
                                     linear_section_index: int,
                                     indent: float=0.1):
        """
        """
        x = self.deformations[linear_section_index]
        k = self.tensions[linear_section_index] / self.deformations[linear_section_index]
        linear_section = lambda x: k * x
        elastic_modulus = (linear_section(x - indent * x) - linear_section(indent * x)) / (x * (1 - 2 * indent))
        # обработать деление на 0
        return elastic_modulus
    
    def _return_elastic_modulus_unload(self,
                                       linear_section_index: int,
                                       transition_index_load: int,
                                       indent: float=0.1):
        """
        """
        x = self.deformations[linear_section_index]
        k = (self.tensions[transition_index_load] - self.tensions[linear_section_index]) / (self.deformations[transition_index_load] - self.deformations[linear_section_index])

        linear_section = lambda x: k * (x - self.deformations[linear_section_index]) + self.tensions[linear_section_index]
        elastic_modulus = (linear_section(self.deformations[transition_index_load] - indent * x) - linear_section(x + indent * x)) / (self.deformations[transition_index_load] - x * (1 + 2 * indent)) 
        # обработать деление на 0
        return elastic_modulus

    def fit(self, 
            admission: float=0.02,
            indent: float=0.1) -> dict:
        """
        """
        index_transition_load = self._find_index_transition_point_load()
        transition_point_load = (self.deformations[index_transition_load], 
                                 self.tensions[index_transition_load])

        index_transition_unload = self._find_index_transition_point_unload()
        transition_point_unload = (self.deformations[index_transition_unload], 
                                   self.tensions[index_transition_unload])

        index_linear_section_load = self._find_index_linear_section_load()
        linear_section_end_point_load = (self.deformations[index_linear_section_load], 
                                         self.tensions[index_linear_section_load])

        index_linear_section_unload = self._find_index_linear_section_unload()
        linear_section_end_point_unload = (self.deformations[index_linear_section_unload], 
                                         self.tensions[index_linear_section_unload])

        limit_strength_load = self._find_limit_strength_load(index_linear_section_load, index_transition_load, admission)
        limit_strength_unload = self._find_limit_strength_unload(index_linear_section_unload, index_transition_load, index_transition_unload, admission)
        
        elastic_modulus_load = self._return_elastic_modulus_load(index_linear_section_load, indent)
        elastic_modulus_unload = self._return_elastic_modulus_unload(index_linear_section_unload, index_transition_load, indent)

        self.info  = {"Transition point load": transition_point_load,
                "Transition point unload": transition_point_unload,
                "End of linear section point load": linear_section_end_point_load,
                "End of linear section point unload": linear_section_end_point_unload,
                "Limit strength point load": limit_strength_load,
                "Limit strength point unload": limit_strength_unload,
                "Elastic modulus load":  elastic_modulus_load,
                "Elastic modulus unload":  elastic_modulus_unload}

        return self.info
    
def make_plot(object: ChartParametrs,
              figsize: Union[list, tuple]=(12, 10),
              linewidth: Union[int, float]=2,
              s: Union[int, float]=8,
              color: str="purple",
              fontsize: int=17) -> plt:

    x = object.deformations
    y = object.tensions
    info = object.fit()

    fig = plt.figure(figsize=figsize)
    plt.plot(x, y, linewidth=linewidth, color=color)
    plt.scatter(x, y, s=s, color=color)
    plt.plot(*info["End of linear section point load"], "o", color="blue")
    plt.plot((x[0], info["End of linear section point load"][0]),
             (y[0], info["End of linear section point load"][1]), linewidth=linewidth, color="blue")

    plt.plot(info["Transition point load"][0], info["Transition point load"][1], "o", color="blue")
    plt.plot(*info["End of linear section point unload"], "o", color="blue")
    plt.plot((info["Transition point load"][0], info["End of linear section point unload"][0]),
             (info["Transition point load"][1], info["End of linear section point unload"][1]), linewidth=linewidth, color="blue")
    
    plt.plot(x[0], y[0], "o", color="red")
    plt.plot(*info["Limit strength point load"], "o", color="red")
    plt.plot((x[0], info["Limit strength point load"][0]),
             (y[0], info["Limit strength point load"][1]), color="red", linewidth=linewidth)

    plt.plot(info["Transition point load"][0], info["Transition point load"][1], "o", color="red")
    plt.plot(*info["Limit strength point unload"], "o", color="red")
    plt.plot((info["Transition point load"][0], info["Limit strength point unload"][0]),
             (info["Transition point load"][1], info["Limit strength point unload"][1]), color="red", linewidth=linewidth)
    
    plt.title(f'E load: {info["Elastic modulus load"]:.2f}, E unload: {info["Elastic modulus load"]:.2f}', fontsize=fontsize)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.grid()
    return plt

if __name__ == "__main__":
    pass