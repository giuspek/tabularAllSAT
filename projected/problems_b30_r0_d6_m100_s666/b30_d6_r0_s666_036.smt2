(set-logic QF_UF)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A12))) (let ((.def_1 (and A6 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A19))) (let ((.def_4 (and .def_3 A16))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not A5))) (let ((.def_7 (or .def_6 A14))) (let ((.def_8 (not .def_7))) (let ((.def_9 (or A18 A17))) (let ((.def_10 (not .def_9))) (let ((.def_11 (or .def_10 .def_8))) (let ((.def_12 (or .def_11 .def_5))) (let ((.def_13 (not .def_12))) (let ((.def_14 (not A3))) (let ((.def_15 (and .def_14 A15))) (let ((.def_16 (not .def_15))) (let ((.def_17 (not A4))) (let ((.def_18 (or A13 .def_17))) (let ((.def_19 (not .def_18))) (let ((.def_20 (or .def_19 .def_16))) (let ((.def_21 (not .def_20))) (let ((.def_22 (and .def_0 .def_6))) (let ((.def_23 (or A24 A18))) (let ((.def_24 (not .def_23))) (let ((.def_25 (and .def_24 .def_22))) (let ((.def_26 (or .def_25 .def_21))) (let ((.def_27 (= .def_26 .def_13))) (let ((.def_28 (not A26))) (let ((.def_29 (or .def_28 A17))) (let ((.def_30 (and A22 A29))) (let ((.def_31 (and .def_30 .def_29))) (let ((.def_32 (not .def_31))) (let ((.def_33 (not A21))) (let ((.def_34 (or .def_33 .def_17))) (let ((.def_35 (not A16))) (let ((.def_36 (and .def_35 .def_0))) (let ((.def_37 (not .def_36))) (let ((.def_38 (and .def_37 .def_34))) (let ((.def_39 (or .def_38 .def_32))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A9))) (let ((.def_42 (not A23))) (let ((.def_43 (and .def_42 .def_41))) (let ((.def_44 (not A14))) (let ((.def_45 (not A13))) (let ((.def_46 (and .def_45 .def_44))) (let ((.def_47 (or .def_46 .def_43))) (let ((.def_48 (not .def_47))) (let ((.def_49 (not A11))) (let ((.def_50 (and A14 .def_49))) (let ((.def_51 (not .def_50))) (let ((.def_52 (or .def_41 A7))) (let ((.def_53 (not .def_52))) (let ((.def_54 (or .def_53 .def_51))) (let ((.def_55 (not .def_54))) (let ((.def_56 (or .def_55 .def_48))) (let ((.def_57 (not .def_56))) (let ((.def_58 (or .def_57 .def_40))) (let ((.def_59 (and .def_58 .def_27))) (let ((.def_60 (not .def_59))) (let ((.def_61 (not A6))) (let ((.def_62 (or .def_61 A1))) (let ((.def_63 (not A28))) (let ((.def_64 (not A1))) (let ((.def_65 (and .def_64 .def_63))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and .def_66 .def_62))) (let ((.def_68 (not A29))) (let ((.def_69 (= A6 .def_68))) (let ((.def_70 (not A10))) (let ((.def_71 (= .def_70 A20))) (let ((.def_72 (or .def_71 .def_69))) (let ((.def_73 (and .def_72 .def_67))) (let ((.def_74 (not .def_73))) (let ((.def_75 (or A5 A8))) (let ((.def_76 (and A20 A21))) (let ((.def_77 (or .def_76 .def_75))) (let ((.def_78 (or A28 .def_41))) (let ((.def_79 (not .def_78))) (let ((.def_80 (and .def_45 A6))) (let ((.def_81 (or .def_80 .def_79))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_82 .def_77))) (let ((.def_84 (= .def_83 .def_74))) (let ((.def_85 (or A12 A23))) (let ((.def_86 (and A11 .def_68))) (let ((.def_87 (or .def_86 .def_85))) (let ((.def_88 (= .def_3 A27))) (let ((.def_89 (not .def_88))) (let ((.def_90 (and A8 .def_70))) (let ((.def_91 (= .def_90 .def_89))) (let ((.def_92 (not .def_91))) (let ((.def_93 (= .def_92 .def_87))) (let ((.def_94 (not .def_93))) (let ((.def_95 (not A2))) (let ((.def_96 (and .def_95 .def_49))) (let ((.def_97 (not A15))) (let ((.def_98 (and .def_97 A17))) (let ((.def_99 (not .def_98))) (let ((.def_100 (and .def_99 .def_96))) (let ((.def_101 (or A22 .def_44))) (let ((.def_102 (and .def_42 .def_45))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_101))) (let ((.def_105 (or .def_104 .def_100))) (let ((.def_106 (or .def_105 .def_94))) (let ((.def_107 (or .def_106 .def_84))) (let ((.def_108 (and .def_107 .def_60))) .def_108))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
