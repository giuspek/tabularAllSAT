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
(assert (let ((.def_0 (and A1 A9))) (let ((.def_1 (not .def_0))) (let ((.def_2 (not A1))) (let ((.def_3 (not A2))) (let ((.def_4 (or .def_3 .def_2))) (let ((.def_5 (not .def_4))) (let ((.def_6 (or .def_5 .def_1))) (let ((.def_7 (not .def_6))) (let ((.def_8 (or A11 A12))) (let ((.def_9 (not A7))) (let ((.def_10 (not A13))) (let ((.def_11 (= .def_10 .def_9))) (let ((.def_12 (= .def_11 .def_8))) (let ((.def_13 (and .def_12 .def_7))) (let ((.def_14 (not A20))) (let ((.def_15 (or .def_14 A5))) (let ((.def_16 (not .def_15))) (let ((.def_17 (and A14 A7))) (let ((.def_18 (not .def_17))) (let ((.def_19 (and .def_18 .def_16))) (let ((.def_20 (= A7 .def_14))) (let ((.def_21 (or A21 A4))) (let ((.def_22 (not .def_21))) (let ((.def_23 (or .def_22 .def_20))) (let ((.def_24 (not .def_23))) (let ((.def_25 (or .def_24 .def_19))) (let ((.def_26 (or .def_25 .def_13))) (let ((.def_27 (not A17))) (let ((.def_28 (or .def_27 A8))) (let ((.def_29 (not .def_28))) (let ((.def_30 (not A6))) (let ((.def_31 (not A22))) (let ((.def_32 (or .def_31 .def_30))) (let ((.def_33 (not .def_32))) (let ((.def_34 (or .def_33 .def_29))) (let ((.def_35 (and .def_31 A22))) (let ((.def_36 (not .def_35))) (let ((.def_37 (not A11))) (let ((.def_38 (not A5))) (let ((.def_39 (and .def_38 .def_37))) (let ((.def_40 (and .def_39 .def_36))) (let ((.def_41 (not .def_40))) (let ((.def_42 (or .def_41 .def_34))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A15))) (let ((.def_45 (not A24))) (let ((.def_46 (or .def_45 .def_44))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or A18 A22))) (let ((.def_49 (not .def_48))) (let ((.def_50 (and .def_49 .def_47))) (let ((.def_51 (or A9 .def_37))) (let ((.def_52 (or A11 .def_3))) (let ((.def_53 (not .def_52))) (let ((.def_54 (or .def_53 .def_51))) (let ((.def_55 (and .def_54 .def_50))) (let ((.def_56 (or .def_55 .def_43))) (let ((.def_57 (or .def_56 .def_26))) (let ((.def_58 (or A18 A20))) (let ((.def_59 (not A19))) (let ((.def_60 (not A14))) (let ((.def_61 (= .def_60 .def_59))) (let ((.def_62 (not .def_61))) (let ((.def_63 (or .def_62 .def_58))) (let ((.def_64 (or .def_31 A19))) (let ((.def_65 (not .def_64))) (let ((.def_66 (and A23 A6))) (let ((.def_67 (and .def_66 .def_65))) (let ((.def_68 (or .def_67 .def_63))) (let ((.def_69 (or A7 A13))) (let ((.def_70 (= .def_45 .def_37))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_71 .def_69))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_38 A20))) (let ((.def_75 (and A16 .def_14))) (let ((.def_76 (and .def_75 .def_74))) (let ((.def_77 (not .def_76))) (let ((.def_78 (or .def_77 .def_73))) (let ((.def_79 (not .def_78))) (let ((.def_80 (or .def_79 .def_68))) (let ((.def_81 (not .def_80))) (let ((.def_82 (not A4))) (let ((.def_83 (not A10))) (let ((.def_84 (and .def_83 .def_82))) (let ((.def_85 (not .def_84))) (let ((.def_86 (not A3))) (let ((.def_87 (and .def_86 A12))) (let ((.def_88 (and .def_87 .def_85))) (let ((.def_89 (and .def_2 .def_37))) (let ((.def_90 (or A15 A16))) (let ((.def_91 (and .def_90 .def_89))) (let ((.def_92 (not .def_91))) (let ((.def_93 (and .def_92 .def_88))) (let ((.def_94 (and .def_45 A13))) (let ((.def_95 (not .def_94))) (let ((.def_96 (or A2 .def_83))) (let ((.def_97 (not .def_96))) (let ((.def_98 (and .def_97 .def_95))) (let ((.def_99 (not A9))) (let ((.def_100 (and .def_44 .def_99))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and .def_2 A2))) (let ((.def_103 (not .def_102))) (let ((.def_104 (and .def_103 .def_101))) (let ((.def_105 (and .def_104 .def_98))) (let ((.def_106 (or .def_105 .def_93))) (let ((.def_107 (or .def_106 .def_81))) (let ((.def_108 (and .def_107 .def_57))) .def_108))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
