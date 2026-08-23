"""Unit tests for the centralized geometryIO surface reader.

These lock in the hardened parsing rules that replaced the previous ad-hoc
parsers, in particular:
  * OBJ ``v `` lines only (``vn``/``vt``/``vp`` excluded).
  * whitespace-robust token splitting (runs of multiple spaces).
  * gzip transparency.
"""
import gzip

import numpy as np
import pytest

import geometryIO


class TestReadObjVertices:
    def test_excludes_vn_vt_and_handles_multispace(self, sample_obj_file, sample_obj_vertices):
        verts = geometryIO.readObjVertices(sample_obj_file)
        assert verts.shape == (4, 3)
        assert np.allclose(verts, np.array(sample_obj_vertices))

    def test_empty_file_returns_empty(self, tmp_path):
        p = tmp_path / "empty.obj"
        p.write_text("# just a comment\n")
        verts = geometryIO.readObjVertices(str(p))
        assert verts.shape == (0, 3)

    def test_gzip_transparent(self, tmp_path, sample_obj_text, sample_obj_vertices):
        p = tmp_path / "sample.obj.gz"
        with gzip.open(p, "wt") as handle:
            handle.write(sample_obj_text)
        verts = geometryIO.readObjVertices(str(p))
        assert np.allclose(verts, np.array(sample_obj_vertices))


class TestReadStlVertices:
    def test_reads_vertices_multispace(self, sample_stl_file, sample_stl_vertices):
        verts = geometryIO.readStlVertices(sample_stl_file)
        assert verts.shape == (3, 3)
        assert np.allclose(verts, np.array(sample_stl_vertices))


class TestReadSurfaceVertices:
    def test_dispatches_obj(self, sample_obj_file, sample_obj_vertices):
        verts = geometryIO.readSurfaceVertices(sample_obj_file)
        assert np.allclose(verts, np.array(sample_obj_vertices))

    def test_dispatches_stl(self, sample_stl_file, sample_stl_vertices):
        verts = geometryIO.readSurfaceVertices(sample_stl_file)
        assert np.allclose(verts, np.array(sample_stl_vertices))


class TestSurfaceBoundingBox:
    def test_obj_bounds(self, sample_obj_file):
        bb = geometryIO.surfaceBoundingBox(sample_obj_file)
        # (minX, minY, minZ, maxX, maxY, maxZ)
        assert bb == pytest.approx((0.0, 0.0, 0.0, 1.0, 1.0, 2.0))

    def test_stl_bounds(self, sample_stl_file):
        bb = geometryIO.surfaceBoundingBox(sample_stl_file)
        assert bb == pytest.approx((0.0, 0.0, 0.0, 1.0, 1.0, 2.0))

    def test_empty_raises(self, tmp_path):
        p = tmp_path / "empty.obj"
        p.write_text("# nothing here\n")
        with pytest.raises(ValueError):
            geometryIO.surfaceBoundingBox(str(p))
