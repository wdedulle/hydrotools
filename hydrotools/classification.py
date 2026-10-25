from osgeo import gdal_array, gdal, ogr, osr
import numpy as np
import shapefile
import pngcanvas

gdal.UseExceptions()


def classify_simple(src, tgt):
    """Classify a remotely sensed image

    Args:
        src (str): Input file name (thermal image)
        tgt (str): Output file name
    """
    # Load the image into numpy using gdal
    srcArr = gdal_array.LoadFile(src)

    # Split the histogram into 20 bins as our classes
    classes = gdal_array.numpy.histogram(srcArr, bins=20)[1]

    # Color look-up table (LUT) - must be len(classes)+1.
    # Specified as R, G, B tuples
    lut = [
        [255, 0, 0],
        [191, 48, 48],
        [166, 0, 0],
        [255, 64, 64],
        [255, 115, 115],
        [255, 116, 0],
        [191, 113, 48],
        [255, 178, 115],
        [0, 153, 153],
        [29, 115, 115],
        [0, 99, 99],
        [166, 75, 0],
        [0, 204, 0],
        [51, 204, 204],
        [255, 150, 64],
        [92, 204, 204],
        [38, 153, 38],
        [0, 133, 0],
        [57, 230, 57],
        [103, 230, 103],
        [184, 138, 0],
    ]

    # Starting value for classification
    start = 1

    # Set up the RGB color JPEG output image
    rgb = gdal_array.numpy.zeros(
        (
            3,
            srcArr.shape[0],
            srcArr.shape[1],
        ),
        gdal_array.numpy.float32,
    )

    # Process all classes and assign colors
    for i in range(len(classes)):
        mask = gdal_array.numpy.logical_and(start <= srcArr, srcArr <= classes[i])
        for j in range(len(lut[i])):
            rgb[j] = gdal_array.numpy.choose(mask, (rgb[j], lut[i][j]))
        start = classes[i] + 1

    # Save the image
    output = gdal_array.SaveArray(
        rgb.astype(gdal_array.numpy.uint8), tgt, format="JPEG"
    )
    output = None


gdal.UseExceptions()

RAMPS = {
    "viridis": [
        (68, 1, 84),
        (59, 82, 139),
        (33, 145, 140),
        (94, 201, 98),
        (253, 231, 37),
    ],
    "blues": [
        (8, 48, 107),
        (33, 113, 181),
        (107, 174, 214),
        (198, 219, 239),
        (247, 251, 255),
    ],
}


def make_lut(n, ramp="viridis"):
    anchors = np.array(RAMPS[ramp], dtype=float)
    pos_anchor = np.linspace(0, 1, len(anchors))
    pos = np.linspace(0, 1, n)
    return (
        np.stack([np.interp(pos, pos_anchor, anchors[:, c]) for c in range(3)], axis=1)
        .round()
        .astype(np.uint8)
    )


def compute_edges(valid, method="quantile", n_classes=10, breaks=None):
    if breaks is not None:
        return np.asarray(breaks, dtype=float)
    if method == "equal":
        return np.linspace(valid.min(), valid.max(), n_classes + 1)
    if method == "quantile":
        return np.unique(np.quantile(valid, np.linspace(0, 1, n_classes + 1)))
    raise ValueError(f"Unknown method: {method}")


def classify_bathymetry(
    src_path,
    tgt_path,
    method="quantile",
    n_classes=10,
    breaks=None,
    ramp="viridis",
    nodata_color=(255, 255, 255),
):
    # 1. Read
    ds = gdal.Open(src_path)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    nodata = band.GetNoDataValue()

    # 2. Nodata-mask (NaN/inf + nodata-values from metadata)
    mask = ~np.isfinite(arr)
    if nodata is not None:
        mask |= arr == nodata
    valid = arr[~mask]

    # 3. Determine class borders and Look Up Table (LUT)  dynamically
    edges = compute_edges(valid, method, n_classes, breaks)
    n = len(edges) - 1
    lut = make_lut(n, ramp)

    # 4. classify: np.digitize gives index 0..n-1
    class_idx = np.clip(np.digitize(arr, edges[1:-1]), 0, n - 1)
    rgb = lut[class_idx]  # (rows, cols, 3) uint8
    rgb[mask] = nodata_color

    rows, cols = arr.shape
    mem = gdal.GetDriverByName("MEM").Create("", cols, rows, 3, gdal.GDT_Byte)
    mem.SetGeoTransform(ds.GetGeoTransform())
    mem.SetProjection(ds.GetProjection())
    for i in range(3):
        mem.GetRasterBand(i + 1).WriteArray(rgb[:, :, i])

    # 6. To JPEG (WORLDFILE=YES writes a .wld besides JPEG)
    gdal.GetDriverByName("JPEG").CreateCopy(
        tgt_path, mem, options=["QUALITY=95", "WORLDFILE=YES"]
    )

    mem = None
    ds = None


def extract_features(
    src: str,
    tgt: str,
    band_index: int = 1,
    nodata_value: int | None = None,
    min_pixels: int = 30,
    layer_name: str = "classes",
):
    """Automatically extract features of a thresholded image to a shapefile

    Args:
        src (str): Thresholded input raster name
        tgt (str): Output shapefile name
    """
    src_ds = gdal.Open(src)
    band = src_ds.GetRasterBand(band_index)
    data = band.ReadAsArray()

    if nodata_value is None:
        nodata_value = band.GetNoDataValue()

    valid = np.ones(data.shape, dtype=bool)
    if nodata_value is not None:
        valid &= data != nodata_value

    # Klasseraster in het geheugen (nodata -> 0, dus klassen moeten >= 1 zijn)
    classes = np.where(valid, data, 0).astype(np.int32)

    mem = gdal.GetDriverByName("MEM")
    cls_ds = mem.Create("", src_ds.RasterXSize, src_ds.RasterYSize, 1, gdal.GDT_Int32)
    cls_ds.SetGeoTransform(src_ds.GetGeoTransform())
    cls_ds.SetProjection(src_ds.GetProjectionRef())
    cls_band = cls_ds.GetRasterBand(1)
    cls_band.WriteArray(classes)

    mask_ds = mem.Create("", src_ds.RasterXSize, src_ds.RasterYSize, 1, gdal.GDT_Byte)
    mask_band = mask_ds.GetRasterBand(1)
    mask_band.WriteArray(valid.astype(np.uint8))

    if min_pixels > 1:
        gdal.SieveFilter(cls_band, mask_band, cls_band, min_pixels, 8)

    driver = ogr.GetDriverByName("ESRI Shapefile")
    shp = driver.CreateDataSource(tgt)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjectionRef())
    layer = shp.CreateLayer(layer_name, srs=srs)
    layer.CreateField(ogr.FieldDefn("CLASS", ogr.OFTInteger))

    gdal.Polygonize(cls_band, mask_band, layer, 0, ["8CONNECTED=8"], None)

    shp = None
    cls_ds = mask_ds = src_ds = None


def rasterize_shapefile_outline(
    shapefile_path: str,
    output_png: str,
    width: int = 800,
    height: int = 600,
):
    # Open the extracted islands
    r = shapefile.Reader(shapefile_path)

    # Setup the world to pixels conversion
    xdist = r.bbox[2] - r.bbox[0]
    ydist = r.bbox[3] - r.bbox[1]

    xratio = width / xdist
    yratio = height / ydist
    polygons = []

    # Loop through all shapes
    for shape in r.shapes():
        # Loop through all parts to catch
        # polygon holes!
        for i in range(len(shape.parts)):
            pixels = []
            pt = None
            if i < len(shape.parts) - 1:
                pt = shape.points[shape.parts[i] : shape.parts[i + 1]]
            else:
                pt = shape.points[shape.parts[i] :]
            for x, y in pt:
                px = int(width - ((r.bbox[2] - x) * xratio))
                py = int((r.bbox[3] - y) * yratio)
                pixels.append([px, py])
            polygons.append(pixels)

    # Set up the output canvas
    c = pngcanvas.PNGCanvas(width, height)

    # Loop through the polygons and draw them
    for p in polygons:
        c.polyline(p)

    # Save the image
    with open(output_png, "wb") as f:
        f.write(c.dump())
        f.close()
