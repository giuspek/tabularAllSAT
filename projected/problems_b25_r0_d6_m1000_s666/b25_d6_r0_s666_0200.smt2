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
(assert (let ((.def_0 (or A13 A3))) (let ((.def_1 (not A16))) (let ((.def_2 (not A19))) (let ((.def_3 (or .def_2 .def_1))) (let ((.def_4 (and .def_3 .def_0))) (let ((.def_5 (not A13))) (let ((.def_6 (= A13 .def_5))) (let ((.def_7 (not .def_6))) (let ((.def_8 (not A8))) (let ((.def_9 (and A7 .def_8))) (let ((.def_10 (and .def_9 .def_7))) (let ((.def_11 (or .def_10 .def_4))) (let ((.def_12 (not A24))) (let ((.def_13 (or .def_12 .def_12))) (let ((.def_14 (and .def_8 A9))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and .def_15 .def_13))) (let ((.def_17 (or A11 .def_12))) (let ((.def_18 (not .def_17))) (let ((.def_19 (or A0 A5))) (let ((.def_20 (not .def_19))) (let ((.def_21 (and .def_20 .def_18))) (let ((.def_22 (and .def_21 .def_16))) (let ((.def_23 (or .def_22 .def_11))) (let ((.def_24 (not .def_23))) (let ((.def_25 (not A21))) (let ((.def_26 (or .def_25 A0))) (let ((.def_27 (not .def_26))) (let ((.def_28 (not A10))) (let ((.def_29 (and A6 .def_28))) (let ((.def_30 (or .def_29 .def_27))) (let ((.def_31 (not .def_30))) (let ((.def_32 (and A10 A15))) (let ((.def_33 (not .def_32))) (let ((.def_34 (and .def_25 A24))) (let ((.def_35 (and .def_34 .def_33))) (let ((.def_36 (or .def_35 .def_31))) (let ((.def_37 (not .def_36))) (let ((.def_38 (not A12))) (let ((.def_39 (and .def_12 .def_38))) (let ((.def_40 (not .def_39))) (let ((.def_41 (or A20 A15))) (let ((.def_42 (not .def_41))) (let ((.def_43 (or .def_42 .def_40))) (let ((.def_44 (not .def_43))) (let ((.def_45 (not A20))) (let ((.def_46 (or .def_25 .def_45))) (let ((.def_47 (not .def_46))) (let ((.def_48 (not A6))) (let ((.def_49 (or A9 .def_48))) (let ((.def_50 (not .def_49))) (let ((.def_51 (and .def_50 .def_47))) (let ((.def_52 (or .def_51 .def_44))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and .def_53 .def_37))) (let ((.def_55 (not .def_54))) (let ((.def_56 (= .def_55 .def_24))) (let ((.def_57 (or A24 A0))) (let ((.def_58 (not .def_57))) (let ((.def_59 (not A15))) (let ((.def_60 (and A19 .def_59))) (let ((.def_61 (or .def_60 .def_58))) (let ((.def_62 (and .def_48 A5))) (let ((.def_63 (not .def_62))) (let ((.def_64 (not A23))) (let ((.def_65 (or A21 .def_64))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and .def_66 .def_63))) (let ((.def_68 (and .def_67 .def_61))) (let ((.def_69 (not A22))) (let ((.def_70 (and .def_69 .def_28))) (let ((.def_71 (and A8 .def_38))) (let ((.def_72 (and .def_71 .def_70))) (let ((.def_73 (not A4))) (let ((.def_74 (and .def_12 .def_73))) (let ((.def_75 (or A3 A5))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_74))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or .def_78 .def_72))) (let ((.def_80 (or .def_79 .def_68))) (let ((.def_81 (not .def_80))) (let ((.def_82 (and A4 A2))) (let ((.def_83 (not .def_82))) (let ((.def_84 (not A0))) (let ((.def_85 (or .def_84 A23))) (let ((.def_86 (and .def_85 .def_83))) (let ((.def_87 (not .def_86))) (let ((.def_88 (not A17))) (let ((.def_89 (or .def_88 .def_64))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or .def_59 A21))) (let ((.def_92 (not .def_91))) (let ((.def_93 (and .def_92 .def_90))) (let ((.def_94 (not .def_93))) (let ((.def_95 (or .def_94 .def_87))) (let ((.def_96 (not .def_95))) (let ((.def_97 (or .def_12 A24))) (let ((.def_98 (not .def_97))) (let ((.def_99 (and A14 A8))) (let ((.def_100 (and .def_99 .def_98))) (let ((.def_101 (not .def_100))) (let ((.def_102 (not A1))) (let ((.def_103 (and .def_102 .def_88))) (let ((.def_104 (or .def_8 A19))) (let ((.def_105 (or .def_104 .def_103))) (let ((.def_106 (not .def_105))) (let ((.def_107 (and .def_106 .def_101))) (let ((.def_108 (or .def_107 .def_96))) (let ((.def_109 (or .def_108 .def_81))) (let ((.def_110 (= .def_109 .def_56))) (let ((.def_111 (not .def_110))) .def_111)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
