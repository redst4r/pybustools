# This is a generated file! Please edit source .ksy file and use kaitai-struct-compiler to rebuild

from pkg_resources import parse_version
import kaitaistruct
from kaitaistruct import KaitaiStruct, KaitaiStream, BytesIO


if parse_version(kaitaistruct.__version__) < parse_version('0.9'):
    raise Exception("Incompatible Kaitai Struct Python API: 0.9 or later is required, but you have %s" % (kaitaistruct.__version__))

class Busfile(KaitaiStruct):
    def __init__(self, _io, _parent=None, _root=None):
        self._io = _io
        self._parent = _parent
        self._root = _root if _root else self
        self._read()

    def _read(self):
        self.header = Busfile.BusHeader(self._io, self, self._root)
        self.records = []
        i = 0
        while not self._io.is_eof():
            self.records.append(Busfile.BusRecord(self._io, self, self._root))
            i += 1


    class BusHeader(KaitaiStruct):
        def __init__(self, _io, _parent=None, _root=None):
            self._io = _io
            self._parent = _parent
            self._root = _root if _root else self
            self._read()

        def _read(self):
            self.magic = self._io.read_bytes(4)
            if not self.magic == b"\x42\x55\x53\x00":
                raise kaitaistruct.ValidationNotEqualError(b"\x42\x55\x53\x00", self.magic, self._io, u"/types/bus_header/seq/0")
            self.version = self._io.read_u4le()
            self.cb_len = self._io.read_u4le()
            self.umi_len = self._io.read_u4le()
            self.freeheader_len = self._io.read_u4le()
            self.freeheader = self._io.read_bytes(self.freeheader_len)


    class BusRecord(KaitaiStruct):
        def __init__(self, _io, _parent=None, _root=None):
            self._io = _io
            self._parent = _parent
            self._root = _root if _root else self
            self._read()

        def _read(self):
            self.cb = self._io.read_u8le()
            self.umi = self._io.read_u8le()
            self.ec = self._io.read_u4le()
            self.count = self._io.read_u4le()
            self.flags = self._io.read_u4le()
            self.nullbytes = self._io.read_bytes(4)
            if not self.nullbytes == b"\x00\x00\x00\x00":
                raise kaitaistruct.ValidationNotEqualError(b"\x00\x00\x00\x00", self.nullbytes, self._io, u"/types/bus_record/seq/5")



