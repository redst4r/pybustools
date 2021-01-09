meta:
  id: busfile
  endian: le
seq:
  - id: header
    type: bus_header
  - id: records
    type: bus_record
    repeat: eos

types:
  bus_header:
    seq:
      - id: magic
        contents: [BUS, 0]
      - id: version
        type: u4
      - id: cb_len
        type: u4
      - id: umi_len
        type: u4
      - id: freeheader_len
        type: u4
      - id: freeheader
        size: freeheader_len
  bus_record:
    seq:
      - id: cb
        type: u8
      - id: umi
        type: u8
      - id: ec
        type: u4
      - id: count
        type: u4
      - id: flags
        type: u4
      - id: nullbytes
        contents: [0, 0, 0, 0]
