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
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (or A2 A15))) (let ((.def_1 (not A10))) (let ((.def_2 (not A7))) (let ((.def_3 (and .def_2 .def_1))) (let ((.def_4 (or .def_3 .def_0))) (let ((.def_5 (not .def_4))) (let ((.def_6 (not A13))) (let ((.def_7 (not A23))) (let ((.def_8 (and .def_7 .def_6))) (let ((.def_9 (not .def_8))) (let ((.def_10 (= A19 A9))) (let ((.def_11 (not .def_10))) (let ((.def_12 (and .def_11 .def_9))) (let ((.def_13 (and .def_12 .def_5))) (let ((.def_14 (not .def_13))) (let ((.def_15 (not A11))) (let ((.def_16 (and .def_15 A19))) (let ((.def_17 (not A3))) (let ((.def_18 (not A5))) (let ((.def_19 (or .def_18 .def_17))) (let ((.def_20 (or .def_19 .def_16))) (let ((.def_21 (not A9))) (let ((.def_22 (and A13 .def_21))) (let ((.def_23 (and A11 A13))) (let ((.def_24 (and .def_23 .def_22))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or .def_25 .def_20))) (let ((.def_27 (or .def_26 .def_14))) (let ((.def_28 (and A4 A17))) (let ((.def_29 (not .def_28))) (let ((.def_30 (not A18))) (let ((.def_31 (or A6 .def_30))) (let ((.def_32 (and .def_31 .def_29))) (let ((.def_33 (not A12))) (let ((.def_34 (and .def_18 .def_33))) (let ((.def_35 (not .def_34))) (let ((.def_36 (and A3 A9))) (let ((.def_37 (not .def_36))) (let ((.def_38 (and .def_37 .def_35))) (let ((.def_39 (or .def_38 .def_32))) (let ((.def_40 (not A21))) (let ((.def_41 (and A11 .def_40))) (let ((.def_42 (not .def_41))) (let ((.def_43 (not A24))) (let ((.def_44 (or A4 .def_43))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_45 .def_42))) (let ((.def_47 (not .def_46))) (let ((.def_48 (= A15 .def_18))) (let ((.def_49 (or A23 A7))) (let ((.def_50 (or .def_49 .def_48))) (let ((.def_51 (and .def_50 .def_47))) (let ((.def_52 (not .def_51))) (let ((.def_53 (or .def_52 .def_39))) (let ((.def_54 (not .def_53))) (let ((.def_55 (or .def_54 .def_27))) (let ((.def_56 (not .def_55))) (let ((.def_57 (not A0))) (let ((.def_58 (and .def_57 A8))) (let ((.def_59 (not .def_58))) (let ((.def_60 (= .def_6 .def_17))) (let ((.def_61 (or .def_60 .def_59))) (let ((.def_62 (not .def_61))) (let ((.def_63 (not A22))) (let ((.def_64 (or A8 .def_63))) (let ((.def_65 (and A20 .def_7))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and .def_66 .def_64))) (let ((.def_68 (or .def_67 .def_62))) (let ((.def_69 (= A11 A2))) (let ((.def_70 (not .def_69))) (let ((.def_71 (and A2 A2))) (let ((.def_72 (and .def_71 .def_70))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_21 .def_18))) (let ((.def_75 (not .def_74))) (let ((.def_76 (or .def_33 A18))) (let ((.def_77 (not .def_76))) (let ((.def_78 (or .def_77 .def_75))) (let ((.def_79 (not .def_78))) (let ((.def_80 (and .def_79 .def_73))) (let ((.def_81 (not .def_80))) (let ((.def_82 (or .def_81 .def_68))) (let ((.def_83 (not .def_82))) (let ((.def_84 (not A20))) (let ((.def_85 (or A4 .def_84))) (let ((.def_86 (not .def_85))) (let ((.def_87 (not A14))) (let ((.def_88 (and .def_87 .def_21))) (let ((.def_89 (not .def_88))) (let ((.def_90 (and .def_89 .def_86))) (let ((.def_91 (not .def_90))) (let ((.def_92 (not A4))) (let ((.def_93 (and .def_92 A5))) (let ((.def_94 (not A16))) (let ((.def_95 (or A19 .def_94))) (let ((.def_96 (not .def_95))) (let ((.def_97 (or .def_96 .def_93))) (let ((.def_98 (not .def_97))) (let ((.def_99 (= .def_98 .def_91))) (let ((.def_100 (not A15))) (let ((.def_101 (= .def_84 .def_100))) (let ((.def_102 (or .def_43 A23))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_101))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and A24 A6))) (let ((.def_107 (or A16 A8))) (let ((.def_108 (not .def_107))) (let ((.def_109 (or .def_108 .def_106))) (let ((.def_110 (or .def_109 .def_105))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_99))) (let ((.def_113 (or .def_112 .def_83))) (let ((.def_114 (or .def_113 .def_56))) (let ((.def_115 (or .def_1 A10))) (let ((.def_116 (not .def_115))) (let ((.def_117 (or A11 A5))) (let ((.def_118 (or .def_117 .def_116))) (let ((.def_119 (not .def_118))) (let ((.def_120 (or A8 .def_33))) (let ((.def_121 (not .def_120))) (let ((.def_122 (not A17))) (let ((.def_123 (and .def_122 .def_122))) (let ((.def_124 (not .def_123))) (let ((.def_125 (and .def_124 .def_121))) (let ((.def_126 (not .def_125))) (let ((.def_127 (and .def_126 .def_119))) (let ((.def_128 (not A2))) (let ((.def_129 (or .def_30 .def_128))) (let ((.def_130 (not .def_129))) (let ((.def_131 (not A8))) (let ((.def_132 (or .def_131 .def_21))) (let ((.def_133 (or .def_132 .def_130))) (let ((.def_134 (= .def_30 A19))) (let ((.def_135 (not .def_134))) (let ((.def_136 (or A23 .def_1))) (let ((.def_137 (not .def_136))) (let ((.def_138 (or .def_137 .def_135))) (let ((.def_139 (and .def_138 .def_133))) (let ((.def_140 (not .def_139))) (let ((.def_141 (or .def_140 .def_127))) (let ((.def_142 (not .def_141))) (let ((.def_143 (or .def_33 .def_84))) (let ((.def_144 (not .def_143))) (let ((.def_145 (and .def_40 A18))) (let ((.def_146 (not .def_145))) (let ((.def_147 (or .def_146 .def_144))) (let ((.def_148 (not .def_147))) (let ((.def_149 (and A2 .def_7))) (let ((.def_150 (and .def_18 .def_100))) (let ((.def_151 (or .def_150 .def_149))) (let ((.def_152 (not .def_151))) (let ((.def_153 (and .def_152 .def_148))) (let ((.def_154 (and A19 A7))) (let ((.def_155 (= .def_122 .def_2))) (let ((.def_156 (and .def_155 .def_154))) (let ((.def_157 (or .def_100 .def_1))) (let ((.def_158 (or .def_128 .def_33))) (let ((.def_159 (not .def_158))) (let ((.def_160 (or .def_159 .def_157))) (let ((.def_161 (not .def_160))) (let ((.def_162 (or .def_161 .def_156))) (let ((.def_163 (not .def_162))) (let ((.def_164 (and .def_163 .def_153))) (let ((.def_165 (not .def_164))) (let ((.def_166 (and .def_165 .def_142))) (let ((.def_167 (and A10 .def_87))) (let ((.def_168 (and .def_2 .def_15))) (let ((.def_169 (not .def_168))) (let ((.def_170 (or .def_169 .def_167))) (let ((.def_171 (not .def_170))) (let ((.def_172 (or A20 .def_84))) (let ((.def_173 (not .def_172))) (let ((.def_174 (or .def_94 .def_94))) (let ((.def_175 (or .def_174 .def_173))) (let ((.def_176 (or .def_175 .def_171))) (let ((.def_177 (or .def_17 A24))) (let ((.def_178 (not .def_177))) (let ((.def_179 (not A1))) (let ((.def_180 (= .def_179 .def_6))) (let ((.def_181 (= .def_180 .def_178))) (let ((.def_182 (and A20 A6))) (let ((.def_183 (not .def_182))) (let ((.def_184 (and A6 .def_1))) (let ((.def_185 (and .def_184 .def_183))) (let ((.def_186 (and .def_185 .def_181))) (let ((.def_187 (or .def_186 .def_176))) (let ((.def_188 (or .def_57 .def_7))) (let ((.def_189 (not A19))) (let ((.def_190 (or .def_189 A14))) (let ((.def_191 (not .def_190))) (let ((.def_192 (and .def_191 .def_188))) (let ((.def_193 (not .def_192))) (let ((.def_194 (or .def_100 .def_84))) (let ((.def_195 (= .def_57 A11))) (let ((.def_196 (= .def_195 .def_194))) (let ((.def_197 (or .def_196 .def_193))) (let ((.def_198 (not .def_197))) (let ((.def_199 (and .def_21 A18))) (let ((.def_200 (not .def_199))) (let ((.def_201 (and .def_21 .def_92))) (let ((.def_202 (not .def_201))) (let ((.def_203 (and .def_202 .def_200))) (let ((.def_204 (and .def_33 .def_100))) (let ((.def_205 (not .def_204))) (let ((.def_206 (or .def_33 .def_30))) (let ((.def_207 (and .def_206 .def_205))) (let ((.def_208 (or .def_207 .def_203))) (let ((.def_209 (not .def_208))) (let ((.def_210 (and .def_209 .def_198))) (let ((.def_211 (not .def_210))) (let ((.def_212 (and .def_211 .def_187))) (let ((.def_213 (or .def_212 .def_166))) (let ((.def_214 (or .def_213 .def_114))) (let ((.def_215 (not .def_214))) (let ((.def_216 (or A2 A3))) (let ((.def_217 (not .def_216))) (let ((.def_218 (and .def_84 A9))) (let ((.def_219 (not .def_218))) (let ((.def_220 (or .def_219 .def_217))) (let ((.def_221 (not .def_220))) (let ((.def_222 (and .def_1 .def_100))) (let ((.def_223 (not .def_222))) (let ((.def_224 (and A1 A12))) (let ((.def_225 (and .def_224 .def_223))) (let ((.def_226 (and .def_225 .def_221))) (let ((.def_227 (and .def_33 .def_15))) (let ((.def_228 (or .def_7 A2))) (let ((.def_229 (= .def_228 .def_227))) (let ((.def_230 (not .def_229))) (let ((.def_231 (and A7 .def_63))) (let ((.def_232 (= .def_30 A17))) (let ((.def_233 (not .def_232))) (let ((.def_234 (and .def_233 .def_231))) (let ((.def_235 (not .def_234))) (let ((.def_236 (or .def_235 .def_230))) (let ((.def_237 (not .def_236))) (let ((.def_238 (or .def_237 .def_226))) (let ((.def_239 (or .def_21 A5))) (let ((.def_240 (not .def_239))) (let ((.def_241 (and .def_84 A23))) (let ((.def_242 (or .def_241 .def_240))) (let ((.def_243 (not .def_242))) (let ((.def_244 (and .def_179 A2))) (let ((.def_245 (not .def_244))) (let ((.def_246 (or .def_122 A1))) (let ((.def_247 (not .def_246))) (let ((.def_248 (and .def_247 .def_245))) (let ((.def_249 (not .def_248))) (let ((.def_250 (= .def_249 .def_243))) (let ((.def_251 (or A19 A14))) (let ((.def_252 (not .def_251))) (let ((.def_253 (or .def_33 .def_189))) (let ((.def_254 (not .def_253))) (let ((.def_255 (or .def_254 .def_252))) (let ((.def_256 (and A1 A23))) (let ((.def_257 (not A6))) (let ((.def_258 (and .def_1 .def_257))) (let ((.def_259 (not .def_258))) (let ((.def_260 (and .def_259 .def_256))) (let ((.def_261 (not .def_260))) (let ((.def_262 (and .def_261 .def_255))) (let ((.def_263 (not .def_262))) (let ((.def_264 (or .def_263 .def_250))) (let ((.def_265 (or .def_264 .def_238))) (let ((.def_266 (and A14 A16))) (let ((.def_267 (or A13 A10))) (let ((.def_268 (not .def_267))) (let ((.def_269 (or .def_268 .def_266))) (let ((.def_270 (not .def_269))) (let ((.def_271 (or A3 .def_189))) (let ((.def_272 (and .def_33 A21))) (let ((.def_273 (not .def_272))) (let ((.def_274 (and .def_273 .def_271))) (let ((.def_275 (not .def_274))) (let ((.def_276 (and .def_275 .def_270))) (let ((.def_277 (not .def_276))) (let ((.def_278 (and .def_2 .def_87))) (let ((.def_279 (or .def_7 .def_43))) (let ((.def_280 (not .def_279))) (let ((.def_281 (or .def_280 .def_278))) (let ((.def_282 (= .def_15 A8))) (let ((.def_283 (or A7 .def_15))) (let ((.def_284 (not .def_283))) (let ((.def_285 (or .def_284 .def_282))) (let ((.def_286 (not .def_285))) (let ((.def_287 (or .def_286 .def_281))) (let ((.def_288 (and .def_287 .def_277))) (let ((.def_289 (not .def_288))) (let ((.def_290 (not .def_227))) (let ((.def_291 (or A1 A3))) (let ((.def_292 (not .def_291))) (let ((.def_293 (and .def_292 .def_290))) (let ((.def_294 (and A1 .def_2))) (let ((.def_295 (= A23 A16))) (let ((.def_296 (or .def_295 .def_294))) (let ((.def_297 (= .def_296 .def_293))) (let ((.def_298 (and A0 A0))) (let ((.def_299 (or A24 .def_189))) (let ((.def_300 (or .def_299 .def_298))) (let ((.def_301 (not .def_300))) (let ((.def_302 (= .def_122 A3))) (let ((.def_303 (not .def_302))) (let ((.def_304 (or .def_131 A5))) (let ((.def_305 (and .def_304 .def_303))) (let ((.def_306 (= .def_305 .def_301))) (let ((.def_307 (not .def_306))) (let ((.def_308 (and .def_307 .def_297))) (let ((.def_309 (or .def_308 .def_289))) (let ((.def_310 (not .def_309))) (let ((.def_311 (and .def_310 .def_265))) (let ((.def_312 (not .def_311))) (let ((.def_313 (or .def_40 .def_179))) (let ((.def_314 (and .def_92 A11))) (let ((.def_315 (or .def_314 .def_313))) (let ((.def_316 (not .def_315))) (let ((.def_317 (= .def_122 .def_1))) (let ((.def_318 (not .def_317))) (let ((.def_319 (= A6 .def_30))) (let ((.def_320 (not .def_319))) (let ((.def_321 (or .def_320 .def_318))) (let ((.def_322 (not .def_321))) (let ((.def_323 (or .def_322 .def_316))) (let ((.def_324 (or A22 .def_131))) (let ((.def_325 (not .def_324))) (let ((.def_326 (and A1 .def_179))) (let ((.def_327 (not .def_326))) (let ((.def_328 (or .def_327 .def_325))) (let ((.def_329 (and .def_18 A14))) (let ((.def_330 (or .def_122 A16))) (let ((.def_331 (and .def_330 .def_329))) (let ((.def_332 (not .def_331))) (let ((.def_333 (and .def_332 .def_328))) (let ((.def_334 (not .def_333))) (let ((.def_335 (= .def_334 .def_323))) (let ((.def_336 (not .def_335))) (let ((.def_337 (or A17 .def_6))) (let ((.def_338 (not .def_337))) (let ((.def_339 (and .def_92 .def_84))) (let ((.def_340 (not .def_339))) (let ((.def_341 (and .def_340 .def_338))) (let ((.def_342 (or A2 .def_6))) (let ((.def_343 (not .def_342))) (let ((.def_344 (or .def_1 A0))) (let ((.def_345 (or .def_344 .def_343))) (let ((.def_346 (or .def_345 .def_341))) (let ((.def_347 (not .def_346))) (let ((.def_348 (= A13 .def_40))) (let ((.def_349 (not .def_348))) (let ((.def_350 (or .def_18 .def_15))) (let ((.def_351 (= .def_350 .def_349))) (let ((.def_352 (not .def_351))) (let ((.def_353 (or .def_87 A15))) (let ((.def_354 (not .def_353))) (let ((.def_355 (and A19 .def_63))) (let ((.def_356 (not .def_355))) (let ((.def_357 (and .def_356 .def_354))) (let ((.def_358 (not .def_357))) (let ((.def_359 (or .def_358 .def_352))) (let ((.def_360 (or .def_359 .def_347))) (let ((.def_361 (not .def_360))) (let ((.def_362 (or .def_361 .def_336))) (let ((.def_363 (not .def_362))) (let ((.def_364 (or A7 .def_40))) (let ((.def_365 (= A9 A7))) (let ((.def_366 (not .def_365))) (let ((.def_367 (or .def_366 .def_364))) (let ((.def_368 (or A2 A1))) (let ((.def_369 (or A23 .def_15))) (let ((.def_370 (not .def_369))) (let ((.def_371 (or .def_370 .def_368))) (let ((.def_372 (not .def_371))) (let ((.def_373 (and .def_372 .def_367))) (let ((.def_374 (or A8 .def_18))) (let ((.def_375 (not .def_374))) (let ((.def_376 (and .def_6 .def_131))) (let ((.def_377 (and .def_376 .def_375))) (let ((.def_378 (or A0 A8))) (let ((.def_379 (not .def_378))) (let ((.def_380 (or A13 A3))) (let ((.def_381 (or .def_380 .def_379))) (let ((.def_382 (or .def_381 .def_377))) (let ((.def_383 (and .def_382 .def_373))) (let ((.def_384 (not .def_383))) (let ((.def_385 (or .def_122 .def_1))) (let ((.def_386 (not .def_385))) (let ((.def_387 (and .def_1 A15))) (let ((.def_388 (= .def_387 .def_386))) (let ((.def_389 (not .def_388))) (let ((.def_390 (or A6 .def_87))) (let ((.def_391 (not .def_390))) (let ((.def_392 (and .def_84 A11))) (let ((.def_393 (not .def_392))) (let ((.def_394 (or .def_393 .def_391))) (let ((.def_395 (or .def_394 .def_389))) (let ((.def_396 (not .def_395))) (let ((.def_397 (and .def_122 A13))) (let ((.def_398 (= .def_17 A23))) (let ((.def_399 (not .def_398))) (let ((.def_400 (= .def_399 .def_397))) (let ((.def_401 (not .def_400))) (let ((.def_402 (and A16 A23))) (let ((.def_403 (and A6 A19))) (let ((.def_404 (and .def_403 .def_402))) (let ((.def_405 (not .def_404))) (let ((.def_406 (= .def_405 .def_401))) (let ((.def_407 (not .def_406))) (let ((.def_408 (or .def_407 .def_396))) (let ((.def_409 (not .def_408))) (let ((.def_410 (and .def_409 .def_384))) (let ((.def_411 (not .def_410))) (let ((.def_412 (and .def_411 .def_363))) (let ((.def_413 (not .def_412))) (let ((.def_414 (and .def_413 .def_312))) (let ((.def_415 (and .def_414 .def_215))) .def_415)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
