(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
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
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A13))) (let ((.def_1 (and .def_0 A23))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A12))) (let ((.def_4 (or .def_3 A16))) (let ((.def_5 (not .def_4))) (let ((.def_6 (or .def_5 .def_2))) (let ((.def_7 (not .def_6))) (let ((.def_8 (not A22))) (let ((.def_9 (not A20))) (let ((.def_10 (or .def_9 .def_8))) (let ((.def_11 (not A9))) (let ((.def_12 (= .def_11 A12))) (let ((.def_13 (and .def_12 .def_10))) (let ((.def_14 (and .def_13 .def_7))) (let ((.def_15 (= A22 A24))) (let ((.def_16 (not .def_15))) (let ((.def_17 (or A15 A20))) (let ((.def_18 (or .def_17 .def_16))) (let ((.def_19 (and A21 A2))) (let ((.def_20 (not .def_19))) (let ((.def_21 (not A17))) (let ((.def_22 (= A21 .def_21))) (let ((.def_23 (not .def_22))) (let ((.def_24 (or .def_23 .def_20))) (let ((.def_25 (not .def_24))) (let ((.def_26 (and .def_25 .def_18))) (let ((.def_27 (not .def_26))) (let ((.def_28 (or .def_27 .def_14))) (let ((.def_29 (not .def_28))) (let ((.def_30 (not A0))) (let ((.def_31 (or A21 .def_30))) (let ((.def_32 (not A5))) (let ((.def_33 (and .def_32 .def_3))) (let ((.def_34 (or .def_33 .def_31))) (let ((.def_35 (not A11))) (let ((.def_36 (and .def_35 A16))) (let ((.def_37 (= .def_30 A1))) (let ((.def_38 (not .def_37))) (let ((.def_39 (or .def_38 .def_36))) (let ((.def_40 (and .def_39 .def_34))) (let ((.def_41 (not A23))) (let ((.def_42 (not A16))) (let ((.def_43 (or .def_42 .def_41))) (let ((.def_44 (not .def_43))) (let ((.def_45 (or .def_41 A9))) (let ((.def_46 (and .def_45 .def_44))) (let ((.def_47 (not A1))) (let ((.def_48 (not A24))) (let ((.def_49 (and .def_48 .def_47))) (let ((.def_50 (not .def_49))) (let ((.def_51 (or .def_42 A14))) (let ((.def_52 (and .def_51 .def_50))) (let ((.def_53 (not .def_52))) (let ((.def_54 (= .def_53 .def_46))) (let ((.def_55 (not .def_54))) (let ((.def_56 (or .def_55 .def_40))) (let ((.def_57 (and .def_56 .def_29))) (let ((.def_58 (= A20 A14))) (let ((.def_59 (or .def_47 .def_9))) (let ((.def_60 (not .def_59))) (let ((.def_61 (or .def_60 .def_58))) (let ((.def_62 (and A2 A7))) (let ((.def_63 (and A21 A8))) (let ((.def_64 (and .def_63 .def_62))) (let ((.def_65 (not .def_64))) (let ((.def_66 (and .def_65 .def_61))) (let ((.def_67 (not .def_66))) (let ((.def_68 (or A4 .def_41))) (let ((.def_69 (not .def_68))) (let ((.def_70 (or .def_48 A23))) (let ((.def_71 (or .def_70 .def_69))) (let ((.def_72 (or A0 .def_8))) (let ((.def_73 (not A4))) (let ((.def_74 (= .def_73 A19))) (let ((.def_75 (and .def_74 .def_72))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_71))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or .def_78 .def_67))) (let ((.def_80 (not A7))) (let ((.def_81 (and .def_8 .def_80))) (let ((.def_82 (not .def_81))) (let ((.def_83 (or .def_80 A20))) (let ((.def_84 (and .def_83 .def_82))) (let ((.def_85 (not A19))) (let ((.def_86 (and A6 .def_85))) (let ((.def_87 (not .def_86))) (let ((.def_88 (not A14))) (let ((.def_89 (or .def_88 A24))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or .def_90 .def_87))) (let ((.def_92 (not .def_91))) (let ((.def_93 (or .def_92 .def_84))) (let ((.def_94 (or .def_85 A1))) (let ((.def_95 (or A8 A10))) (let ((.def_96 (not .def_95))) (let ((.def_97 (and .def_96 .def_94))) (let ((.def_98 (or A12 A15))) (let ((.def_99 (and A15 A2))) (let ((.def_100 (not .def_99))) (let ((.def_101 (or .def_100 .def_98))) (let ((.def_102 (and .def_101 .def_97))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_93))) (let ((.def_105 (= .def_104 .def_79))) (let ((.def_106 (not .def_105))) (let ((.def_107 (or .def_106 .def_57))) (let ((.def_108 (not .def_107))) .def_108))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
