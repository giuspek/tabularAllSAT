(set-logic QF_UF)
(declare-fun A0 () Bool)
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
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (or A22 A8))) (let ((.def_1 (not .def_0))) (let ((.def_2 (not A20))) (let ((.def_3 (and .def_2 A11))) (let ((.def_4 (and .def_3 .def_1))) (let ((.def_5 (not .def_4))) (let ((.def_6 (not A24))) (let ((.def_7 (or A25 .def_6))) (let ((.def_8 (not .def_7))) (let ((.def_9 (and A19 A0))) (let ((.def_10 (and .def_9 .def_8))) (let ((.def_11 (or .def_10 .def_5))) (let ((.def_12 (or A6 .def_6))) (let ((.def_13 (not .def_12))) (let ((.def_14 (not A3))) (let ((.def_15 (or .def_14 A13))) (let ((.def_16 (not .def_15))) (let ((.def_17 (or .def_16 .def_13))) (let ((.def_18 (not .def_17))) (let ((.def_19 (not A8))) (let ((.def_20 (not A9))) (let ((.def_21 (= .def_20 .def_19))) (let ((.def_22 (not .def_21))) (let ((.def_23 (not A4))) (let ((.def_24 (= .def_23 A16))) (let ((.def_25 (not .def_24))) (let ((.def_26 (and .def_25 .def_22))) (let ((.def_27 (and .def_26 .def_18))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_11))) (let ((.def_30 (and A2 A5))) (let ((.def_31 (not .def_30))) (let ((.def_32 (or A12 .def_19))) (let ((.def_33 (not .def_32))) (let ((.def_34 (and .def_33 .def_31))) (let ((.def_35 (not A1))) (let ((.def_36 (or A23 .def_35))) (let ((.def_37 (not A21))) (let ((.def_38 (and A8 .def_37))) (let ((.def_39 (or .def_38 .def_36))) (let ((.def_40 (or .def_39 .def_34))) (let ((.def_41 (or .def_35 A19))) (let ((.def_42 (not .def_41))) (let ((.def_43 (not A19))) (let ((.def_44 (or A13 .def_43))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_45 .def_42))) (let ((.def_47 (not .def_46))) (let ((.def_48 (not A11))) (let ((.def_49 (or .def_48 A17))) (let ((.def_50 (not .def_49))) (let ((.def_51 (not A7))) (let ((.def_52 (and .def_51 A13))) (let ((.def_53 (or .def_52 .def_50))) (let ((.def_54 (or .def_53 .def_47))) (let ((.def_55 (not .def_54))) (let ((.def_56 (and .def_55 .def_40))) (let ((.def_57 (not .def_56))) (let ((.def_58 (and .def_57 .def_29))) (let ((.def_59 (and .def_37 .def_19))) (let ((.def_60 (not .def_59))) (let ((.def_61 (= A20 A22))) (let ((.def_62 (not .def_61))) (let ((.def_63 (or .def_62 .def_60))) (let ((.def_64 (= .def_14 .def_35))) (let ((.def_65 (not .def_64))) (let ((.def_66 (not A27))) (let ((.def_67 (and A1 .def_66))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_68 .def_65))) (let ((.def_70 (not .def_69))) (let ((.def_71 (or .def_70 .def_63))) (let ((.def_72 (or A23 .def_66))) (let ((.def_73 (= A12 .def_6))) (let ((.def_74 (= .def_73 .def_72))) (let ((.def_75 (or A3 .def_35))) (let ((.def_76 (not .def_75))) (let ((.def_77 (not A26))) (let ((.def_78 (or .def_77 .def_37))) (let ((.def_79 (or .def_78 .def_76))) (let ((.def_80 (or .def_79 .def_74))) (let ((.def_81 (and .def_80 .def_71))) (let ((.def_82 (not A12))) (let ((.def_83 (and .def_82 .def_43))) (let ((.def_84 (and A28 A19))) (let ((.def_85 (and .def_84 .def_83))) (let ((.def_86 (not A0))) (let ((.def_87 (and .def_86 A19))) (let ((.def_88 (not A28))) (let ((.def_89 (or A29 .def_88))) (let ((.def_90 (not .def_89))) (let ((.def_91 (and .def_90 .def_87))) (let ((.def_92 (and .def_91 .def_85))) (let ((.def_93 (and A21 .def_43))) (let ((.def_94 (not .def_93))) (let ((.def_95 (not A10))) (let ((.def_96 (and A10 .def_95))) (let ((.def_97 (and .def_96 .def_94))) (let ((.def_98 (not .def_97))) (let ((.def_99 (or A23 A17))) (let ((.def_100 (not .def_99))) (let ((.def_101 (or A28 .def_86))) (let ((.def_102 (or .def_101 .def_100))) (let ((.def_103 (or .def_102 .def_98))) (let ((.def_104 (not .def_103))) (let ((.def_105 (and .def_104 .def_92))) (let ((.def_106 (not .def_105))) (let ((.def_107 (or .def_106 .def_81))) (let ((.def_108 (not .def_107))) (let ((.def_109 (or .def_108 .def_58))) (let ((.def_110 (not .def_109))) .def_110))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)