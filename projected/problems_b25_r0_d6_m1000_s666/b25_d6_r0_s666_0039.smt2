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
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
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
(assert (let ((.def_0 (not A0))) (let ((.def_1 (or A16 .def_0))) (let ((.def_2 (not A11))) (let ((.def_3 (and .def_2 A20))) (let ((.def_4 (= .def_3 .def_1))) (let ((.def_5 (not .def_4))) (let ((.def_6 (and A20 A20))) (let ((.def_7 (not .def_6))) (let ((.def_8 (or A3 A20))) (let ((.def_9 (not .def_8))) (let ((.def_10 (or .def_9 .def_7))) (let ((.def_11 (not .def_10))) (let ((.def_12 (and .def_11 .def_5))) (let ((.def_13 (not .def_12))) (let ((.def_14 (or A21 A8))) (let ((.def_15 (not A19))) (let ((.def_16 (and A16 .def_15))) (let ((.def_17 (or .def_16 .def_14))) (let ((.def_18 (not .def_17))) (let ((.def_19 (not A5))) (let ((.def_20 (or .def_2 .def_19))) (let ((.def_21 (not .def_20))) (let ((.def_22 (and A20 A10))) (let ((.def_23 (or .def_22 .def_21))) (let ((.def_24 (and .def_23 .def_18))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or .def_25 .def_13))) (let ((.def_27 (and .def_2 A19))) (let ((.def_28 (not .def_27))) (let ((.def_29 (not A10))) (let ((.def_30 (or A8 .def_29))) (let ((.def_31 (or .def_30 .def_28))) (let ((.def_32 (not A12))) (let ((.def_33 (or .def_32 A13))) (let ((.def_34 (not .def_33))) (let ((.def_35 (or .def_32 .def_15))) (let ((.def_36 (not .def_35))) (let ((.def_37 (= .def_36 .def_34))) (let ((.def_38 (or .def_37 .def_31))) (let ((.def_39 (not .def_38))) (let ((.def_40 (not A16))) (let ((.def_41 (or .def_40 .def_19))) (let ((.def_42 (not .def_41))) (let ((.def_43 (or A13 A0))) (let ((.def_44 (and .def_43 .def_42))) (let ((.def_45 (or A18 A23))) (let ((.def_46 (not .def_45))) (let ((.def_47 (not A24))) (let ((.def_48 (not A17))) (let ((.def_49 (and .def_48 .def_47))) (let ((.def_50 (not .def_49))) (let ((.def_51 (and .def_50 .def_46))) (let ((.def_52 (and .def_51 .def_44))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and .def_53 .def_39))) (let ((.def_55 (not .def_54))) (let ((.def_56 (or .def_55 .def_26))) (let ((.def_57 (not .def_56))) (let ((.def_58 (not A2))) (let ((.def_59 (and A6 .def_58))) (let ((.def_60 (not .def_59))) (let ((.def_61 (or A18 A4))) (let ((.def_62 (and .def_61 .def_60))) (let ((.def_63 (or A19 A10))) (let ((.def_64 (not A22))) (let ((.def_65 (= A13 .def_64))) (let ((.def_66 (or .def_65 .def_63))) (let ((.def_67 (or .def_66 .def_62))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_40 .def_58))) (let ((.def_70 (not .def_69))) (let ((.def_71 (not A6))) (let ((.def_72 (= .def_19 .def_71))) (let ((.def_73 (not .def_72))) (let ((.def_74 (or .def_73 .def_70))) (let ((.def_75 (not A20))) (let ((.def_76 (or A1 .def_75))) (let ((.def_77 (not A23))) (let ((.def_78 (or .def_77 .def_77))) (let ((.def_79 (not .def_78))) (let ((.def_80 (and .def_79 .def_76))) (let ((.def_81 (not .def_80))) (let ((.def_82 (or .def_81 .def_74))) (let ((.def_83 (not .def_82))) (let ((.def_84 (or .def_83 .def_68))) (let ((.def_85 (and A10 A3))) (let ((.def_86 (or A11 A2))) (let ((.def_87 (or .def_86 .def_85))) (let ((.def_88 (or .def_15 A3))) (let ((.def_89 (or .def_32 A3))) (let ((.def_90 (not .def_89))) (let ((.def_91 (= .def_90 .def_88))) (let ((.def_92 (or .def_91 .def_87))) (let ((.def_93 (not A21))) (let ((.def_94 (or .def_40 .def_93))) (let ((.def_95 (not .def_94))) (let ((.def_96 (and A15 A10))) (let ((.def_97 (or .def_96 .def_95))) (let ((.def_98 (not .def_97))) (let ((.def_99 (and A24 A23))) (let ((.def_100 (not .def_99))) (let ((.def_101 (or A7 .def_77))) (let ((.def_102 (and .def_101 .def_100))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_98))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_92))) (let ((.def_107 (not .def_106))) (let ((.def_108 (= .def_107 .def_84))) (let ((.def_109 (and .def_108 .def_57))) .def_109)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
