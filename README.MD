# Модуль RayTracing
## Описание
Этот модуль предоставляет класс `RayTracing`, предназначенный для трассировки
лучей в сейсмических моделях. Класс позволяет создавать графики лучей,
проходящих через слои с различными параметрами.

## Использование
### 1. Установка зависимостей
Установка зависимостей выполняется командой:
```commandline
poetry install
```

### 2. Входные данные
`incid_angle` - это изначальный угол падения в градусах
`source_depth` - это альтитуда источника луча в метрах
`seismic_model` - это сейсмические данные: альтитуда подошвы пласта в
метрах, альтитуда кровли в метрах, скорость луча в среде в м/с

### 3. Формат входных данных
`incid_angle` - `int` или float
`source_depth` - int или float
`seismic_model` - список списков, в котором значения начинаются с самого
глубинного, то есть с подошвы самого нижнего пласта. Списки внутри
основного списка содержат три элемента в следующей последовательности:
[альтитуда подошвы, альтитуда кровли, скорость луча в среде]

### 4. Запуск модуля
Создается экземпляр класса `RayTracing` c входными данными:
```commandline
ray_tracing = Rayracing
    seismic_model=[
        [-2500, -2200, 3000],
        [-2200, -1800, 2500],
        [-1800, -1000, 2000]
    ],
    incid_angle=39.85,
    source_depth=-2300
```
Вызывается метод `trace_ray`:
```commandline
ray_tracing.trace_ray()
```
Запускается интерпретатор