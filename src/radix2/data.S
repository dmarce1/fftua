#define       M_SQRT1_2      0.70710678118654752440


              .global        tw45
              .global        none
              .global        two
              .global        laneswap16
              .global        laneswap32
              .global        cos45
              .global        cos135
              .global        sin45
              .global        sin135
              .global        one
              .global        half
              .global        zero

              .data

              .align         32
tw45:         .double        M_SQRT1_2
              .double        M_SQRT1_2
              .double        M_SQRT1_2
              .double        M_SQRT1_2
none:         .double        -1.0
              .double        -1.0
              .double        -1.0
              .double        -1.0
one:          .double        1.0
              .double        1.0
              .double        1.0
              .double        1.0
half:         .double        0.5
              .double        0.5
              .double        0.5
              .double        0.5
two:          .double        2.0
              .double        2.0
              .double        2.0
              .double        2.0
zero:         .double        0.0
              .double        0.0
              .double        0.0
              .double        0.0
cos45:        .double        M_SQRT1_2
              .double        M_SQRT1_2
              .double        M_SQRT1_2
              .double        M_SQRT1_2
cos135:       .double        -M_SQRT1_2
              .double        -M_SQRT1_2
              .double        -M_SQRT1_2
              .double        -M_SQRT1_2
sin45:        .double        -M_SQRT1_2
              .double        -M_SQRT1_2
              .double        -M_SQRT1_2
              .double        -M_SQRT1_2
sin135:       .double        -M_SQRT1_2
              .double        -M_SQRT1_2
              .double        -M_SQRT1_2
              .double        -M_SQRT1_2

              .align         16
laneswap16:   .byte          8
              .byte          9
              .byte          10
              .byte          11
              .byte          12
              .byte          13
              .byte          14
              .byte          15
              .byte          0
              .byte          1
              .byte          2
              .byte          3
              .byte          4
              .byte          5
              .byte          6
              .byte          7


