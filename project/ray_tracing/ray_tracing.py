"""Module with class for ray tracing."""

from math import fabs, asin, degrees, radians, sin, tan, sqrt
from typing import List, Union

from matplotlib import pyplot


class RayTracing:
    def __init__(self, seismic_model: List[list],
                 incid_angle: float, source_depth: float):
        """Initializing method.

        Args:
            seismic_model: seismic model
            incid_angle: incidance angle
            source_depth: source depth
        """
        self.seismic_model = seismic_model
        self.incid_angle = incid_angle
        self.source_depth = source_depth

    def check_input_data_correctness(self) -> Union[callable, str]:
        """Check input data correctness.

        Returns: trace_ray method or message
        """
        if not self.is_layer_data_count_correct():
            return 'Layer data is incorrect'
        if not self.is_input_data_integer_or_float():
            return 'Input data is not integer or float'
        if not self.is_speed_data_correct():
            return 'Speed value is not correct'
        if not self.is_source_occur_in_seismic_model():
            return "Source dont occur in seismic model."
        return self.trace_ray()

    def trace_ray(self) -> callable:
        """Method for ray tracing.

        Returns: build_graph method
        """
        layers_values, label_values = self.calculate_ray_tracing_values()
        coords = self.create_coordinates(layers_values)
        self.build_graph(coords)
        self.create_labels(coords=coords, label_values=label_values)

    def calculate_ray_tracing_values(self) -> List[list]:
        """Calculate ray tracing values (ray x-length, ray y-length,
        incidance angle, ray time).

        Returns: list of ray tracing values
        """

        layers_values = []
        label_values = []
        cur_angle = 0
        existing_counter = 0

        for i in range(len(self.seismic_model)):
            is_exist = self.is_layer_exist_source(layer_number=i)

            thickness = self.get_thickness(
                layer_bottom=self.seismic_model[i][1],
                layer_top=self.seismic_model[i][0]
            )
            if is_exist and existing_counter == 0:
                cur_angle = self.incid_angle
                thickness = self.get_thickness(
                    layer_bottom=self.seismic_model[i][1],
                    layer_top=self.source_depth
                )
                existing_counter += 1
            if not is_exist and existing_counter == 0:
                continue

            x_length = self.find_x_length(
                    triangle_angle=cur_angle, thickness=thickness
                )
            ray_length = self.find_ray_length(
                thickness=thickness,
                x_length=x_length
            )
            time = self.count_ray_duration(
                ray_length=ray_length,
                layer_speed=self.seismic_model[i][2]
            )
            layers_values.append([x_length, thickness])
            label_values.append([cur_angle, time])

            if i == len(self.seismic_model) - 1:
                continue

            cur_angle = self.count_angle(
                cur_angle=cur_angle,
                lower_layer_speed=self.seismic_model[i][2],
                upper_layer_speed=self.seismic_model[i+1][2]
            )
            if cur_angle == -2:
                break
        return [layers_values, label_values]

    def is_input_data_integer_or_float(self) -> bool:
        """Check input data type.

        Returns: True if correct, otherwise - False
        """
        for layer in self.seismic_model:
            for i in range(3):
                if type(layer[i]) is int or float:
                    continue
                else:
                    return False
        if type(self.incid_angle) is not int or float:
            return False
        if type(self.source_depth) is not int or float:
            return False
        return True

    def is_speed_data_correct(self) -> bool:
        """Check speed data correctness.

        Returns: True if correct, otherwise - False
        """
        for layer in self.seismic_model:
            if layer[2] > 0:
                continue
            else:
                return False
        return True

    def is_layer_data_count_correct(self) -> bool:
        """Check layer data full (3 value).

        Returns: True if correct, otherwise - False
        """
        for layer in self.seismic_model:
            if len(layer) == 3:
                continue
            else:
                return False
        return True

    def is_source_occur_in_seismic_model(self) -> bool:
        """Check source existing in seismic model.

        Returns: source existing
        """
        for layer in self.seismic_model:
            if self.source_depth in range(layer[0], layer[1]):
                return True
        return False
 
    def get_thickness(self, layer_bottom: float, layer_top: float) -> float:
        """Count thickness of layer.

        Args:
            layer_bottom: layer bottom altitude
            layer_top: layer top altitude

        Returns: layer thickness
        """
        return fabs(layer_bottom - layer_top)

    def is_layer_exist_source(self, layer_number: int) -> bool:
        """Check that source in current layer.

        Args:
            layer_number: layer number

        Returns: True if source in current layer, otherwise - False
        """
        if self.source_depth not in range(
            self.seismic_model[layer_number][0],
            self.seismic_model[layer_number][1]
        ):
            return False
        return True
 
    def count_angle(self, cur_angle: float, lower_layer_speed: float,
                    upper_layer_speed: float) -> float:
        """Count angle of ray in second environment.

        Args:
            cur_angle: current angle
            lower_layer_speed: speed of ray in first layer
            upper_layer_speed: speed of ray in second layer

        Returns: angle of refraction
        """
        sin_beta = (
            sin(radians(cur_angle)) * (upper_layer_speed / lower_layer_speed)
        )
        if fabs(sin_beta) > 1:
            return -2
        return degrees(asin(sin_beta))
 
    def find_x_length(self, triangle_angle: float, thickness: float) -> float:
        """Finds the value of the x_length of a right triangle consisting
        of thickness and ray.

        Args:
            triangle_angle: layer triangle angle
            thickness: thickness of layer

        Returns: value of x_length
        """
        x_length = tan(radians(triangle_angle)) * thickness
        return x_length
 
    def find_ray_length(self, thickness: float, x_length: float) -> float:
        """Find triangle ray_length.

        Args:
            thickness: layer thickness(a-x_length)
            x_length: b-x_length
        """
        ray_length = thickness ** 2 + x_length ** 2
        return sqrt(ray_length)
 
    def count_ray_duration(self, ray_length: float, layer_speed: float
                           ) -> float:
        """Count ray duration.

        Args:
            ray_length: triangle ray_length
            layer_speed: ray speed in layer

        Returns: time in second
        """
        time = ray_length / layer_speed
        return time

    def create_coordinates(self, layers_values: List[list]) -> list:
        """Create coordinates from x_length and thickness.

        Args:
            layers_values: x_length and thickness

        Returns: list with coordinates tuples
        """
        x = [0]
        y = [self.source_depth]
        for i in range(len(layers_values)):
            x_length = layers_values[i][0] + x[i]
            y_length = layers_values[i][1] + y[i]
            x.append(x_length)
            y.append(y_length)
        coords = [x, y]
        return coords
 
    def build_graph(self, coords: List[list]) -> callable:
        """Build graph from coordinates.

        Args:
            coords: coordinates

        Returns: create_label method
        """
        for y_axis in coords[1]:
            pyplot.axhline(y=y_axis, color='r', linestyle='--')
        pyplot.plot(coords[0], coords[1])
        pyplot.scatter(coords[0], coords[1])
 
    def create_labels(self, coords: List[list], label_values: List[list]):
        """Create labels on graph (x-coordinate, y-coordinate,
        incidance angle, ray time).
        """
        for i in range(len(coords[0])):
            if i != len(coords[0]) - 1:
                label = (f'x:{coords[0][i]}\ny:{coords[1][i]}\n'
                         f'incidance angle:{label_values[i][0]}\n'
                         f'time:{label_values[i][1]}')
            else:
                label = f'x:{coords[0][i]}\ny:{coords[1][i]}'
            pyplot.annotate(text=label, xy=(coords[0][i], coords[1][i]))

        pyplot.show()


seismic_model = [
    (-2500, -2200, 3000),
    (-2200, -1800, 2500),
    (-1800, -1000, 2000),
    (-1000, -500, 1500),
    (-500, 200, 900),
    (200, 500, 600)
]
# seismic_model = [(-12, -9, 2), (-9, -6, 3), (-6, -3, 4), (-3, 0, 5)]
q = RayTracing(
    seismic_model=seismic_model,
    incid_angle=39.85,
    source_depth=-2300
)
q.trace_ray()