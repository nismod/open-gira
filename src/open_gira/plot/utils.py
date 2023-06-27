"""
Plotting utilities.
"""


def figure_size(
    min_x: float,
    min_y: float,
    max_x: float,
    max_y: float,
    max_plot_width_in: float = 16,
    max_plot_height_in: float = 9
) -> tuple[float, float]:
    """
    Given bounding box, calculate size of figure in inches for equal aspect ratio.
    Will not exceed `max_plot_width_in` or `max_plot_height_in`.

    Example usage:
    When used in conjunction with a GeoDataFrame named `gdf`:
    ```
    f, ax = plt.subplots(figsize=figure_size(*gdf.total_bounds))
    ```

    Args:
        min_x: Minimum x value in data coordinates
        min_y: Minimum y value in data coordinates
        max_x: Maximum x value in data coordinates
        max_y: Maximum y value in data coordinates

    Returns:
        width in inches, height in inches
    """
    x_span = max_x - min_x
    y_span = max_y - min_y
    aspect_ratio = y_span / x_span
    max_plot_width_in = 16
    max_plot_height_in = 9

    if max_plot_width_in * aspect_ratio < max_plot_height_in:
        # tall
        x_in = max_plot_width_in
        y_in = max_plot_width_in * aspect_ratio

    else:
        # wide
        x_in = max_plot_height_in / aspect_ratio
        y_in = max_plot_height_in

    return x_in, y_in
