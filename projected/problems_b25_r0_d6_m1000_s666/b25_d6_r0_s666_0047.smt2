(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
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
(assert (let ((.def_0 (not A18))) (let ((.def_1 (and A1 .def_0))) (let ((.def_2 (not A14))) (let ((.def_3 (or .def_2 A4))) (let ((.def_4 (or .def_3 .def_1))) (let ((.def_5 (not A19))) (let ((.def_6 (or A3 .def_5))) (let ((.def_7 (= A0 A20))) (let ((.def_8 (not .def_7))) (let ((.def_9 (and .def_8 .def_6))) (let ((.def_10 (not .def_9))) (let ((.def_11 (and .def_10 .def_4))) (let ((.def_12 (not .def_11))) (let ((.def_13 (not A4))) (let ((.def_14 (not A7))) (let ((.def_15 (or .def_14 .def_13))) (let ((.def_16 (= A10 A5))) (let ((.def_17 (or .def_16 .def_15))) (let ((.def_18 (not A5))) (let ((.def_19 (not A23))) (let ((.def_20 (and .def_19 .def_18))) (let ((.def_21 (not .def_20))) (let ((.def_22 (not A2))) (let ((.def_23 (and .def_22 A11))) (let ((.def_24 (not .def_23))) (let ((.def_25 (or .def_24 .def_21))) (let ((.def_26 (and .def_25 .def_17))) (let ((.def_27 (and .def_26 .def_12))) (let ((.def_28 (not .def_27))) (let ((.def_29 (= A14 .def_2))) (let ((.def_30 (not .def_29))) (let ((.def_31 (not A17))) (let ((.def_32 (not A13))) (let ((.def_33 (or .def_32 .def_31))) (let ((.def_34 (or .def_33 .def_30))) (let ((.def_35 (not A11))) (let ((.def_36 (or A11 .def_35))) (let ((.def_37 (not .def_36))) (let ((.def_38 (and A14 A6))) (let ((.def_39 (and .def_38 .def_37))) (let ((.def_40 (or .def_39 .def_34))) (let ((.def_41 (not .def_40))) (let ((.def_42 (not A16))) (let ((.def_43 (and .def_42 .def_2))) (let ((.def_44 (not .def_43))) (let ((.def_45 (or A11 A18))) (let ((.def_46 (not .def_45))) (let ((.def_47 (or .def_46 .def_44))) (let ((.def_48 (not .def_47))) (let ((.def_49 (and A14 .def_5))) (let ((.def_50 (not .def_49))) (let ((.def_51 (not A3))) (let ((.def_52 (and A22 .def_51))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and .def_53 .def_50))) (let ((.def_55 (and .def_54 .def_48))) (let ((.def_56 (and .def_55 .def_41))) (let ((.def_57 (or .def_56 .def_28))) (let ((.def_58 (not .def_57))) (let ((.def_59 (not A20))) (let ((.def_60 (or A24 .def_59))) (let ((.def_61 (not .def_60))) (let ((.def_62 (and A7 A9))) (let ((.def_63 (not .def_62))) (let ((.def_64 (and .def_63 .def_61))) (let ((.def_65 (not .def_64))) (let ((.def_66 (= .def_14 .def_22))) (let ((.def_67 (not .def_66))) (let ((.def_68 (not A6))) (let ((.def_69 (and .def_68 A9))) (let ((.def_70 (not .def_69))) (let ((.def_71 (and .def_70 .def_67))) (let ((.def_72 (not .def_71))) (let ((.def_73 (and .def_72 .def_65))) (let ((.def_74 (not .def_73))) (let ((.def_75 (= A22 .def_32))) (let ((.def_76 (not .def_75))) (let ((.def_77 (not A21))) (let ((.def_78 (or .def_77 .def_68))) (let ((.def_79 (not .def_78))) (let ((.def_80 (= .def_79 .def_76))) (let ((.def_81 (not .def_80))) (let ((.def_82 (not A24))) (let ((.def_83 (and A7 .def_82))) (let ((.def_84 (not .def_83))) (let ((.def_85 (and A23 A19))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and .def_86 .def_84))) (let ((.def_88 (and .def_87 .def_81))) (let ((.def_89 (not .def_88))) (let ((.def_90 (or .def_89 .def_74))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and .def_32 .def_19))) (let ((.def_93 (not .def_92))) (let ((.def_94 (or .def_35 A0))) (let ((.def_95 (not .def_94))) (let ((.def_96 (or .def_95 .def_93))) (let ((.def_97 (not A9))) (let ((.def_98 (or .def_82 .def_97))) (let ((.def_99 (and .def_18 A1))) (let ((.def_100 (not .def_99))) (let ((.def_101 (= .def_100 .def_98))) (let ((.def_102 (or .def_101 .def_96))) (let ((.def_103 (= A15 A3))) (let ((.def_104 (not A10))) (let ((.def_105 (or A23 .def_104))) (let ((.def_106 (not .def_105))) (let ((.def_107 (and .def_106 .def_103))) (let ((.def_108 (not .def_107))) (let ((.def_109 (not A0))) (let ((.def_110 (and .def_18 .def_109))) (let ((.def_111 (and A3 .def_22))) (let ((.def_112 (or .def_111 .def_110))) (let ((.def_113 (not .def_112))) (let ((.def_114 (and .def_113 .def_108))) (let ((.def_115 (not .def_114))) (let ((.def_116 (or .def_115 .def_102))) (let ((.def_117 (not .def_116))) (let ((.def_118 (or .def_117 .def_91))) (let ((.def_119 (or .def_118 .def_58))) .def_119)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
