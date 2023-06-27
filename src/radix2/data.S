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
              .global        UR0
              .global        UR1
              .global        UR2
              .global        UR3
              .global        UI0
              .global        UI1
              .global        UI2
              .global        UI3
              .global        COS1
              .global        SIN1
              .global        COS2
              .global        SIN2
              .global        COS3
              .global        SIN3
              .global        COS4
              .global        SIN4
              .global        COS5
              .global        SIN5
              .global        COS6
              .global        SIN6
              .global        COS7
              .global        SIN7

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
UR0:          .double        0.0, 0.0, 0.0, 0.0
UR1:          .double        0.0, 0.0, 0.0, 0.0
UR2:          .double        0.0, 0.0, 0.0, 0.0
UR3:          .double        0.0, 0.0, 0.0, 0.0
UI0:          .double        0.0, 0.0, 0.0, 0.0
UI1:          .double        0.0, 0.0, 0.0, 0.0
UI2:          .double        0.0, 0.0, 0.0, 0.0
UI3:          .double        0.0, 0.0, 0.0, 0.0
COS1:         .double        0.0, 0.0, 0.0, 0.0
SIN1:         .double        0.0, 0.0, 0.0, 0.0
COS2:         .double        0.0, 0.0, 0.0, 0.0
SIN2:         .double        0.0, 0.0, 0.0, 0.0
COS3:         .double        0.0, 0.0, 0.0, 0.0
SIN3:         .double        0.0, 0.0, 0.0, 0.0
COS4:         .double        0.0, 0.0, 0.0, 0.0
SIN4:         .double        0.0, 0.0, 0.0, 0.0
COS5:         .double        0.0, 0.0, 0.0, 0.0
SIN5:         .double        0.0, 0.0, 0.0, 0.0
COS6:         .double        0.0, 0.0, 0.0, 0.0
SIN6:         .double        0.0, 0.0, 0.0, 0.0
COS7:         .double        0.0, 0.0, 0.0, 0.0
SIN7:         .double        0.0, 0.0, 0.0, 0.0


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


